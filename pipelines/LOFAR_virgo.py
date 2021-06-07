#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# perform self-calibration on a group of SBs concatenated in TCs.
# pipeline for a bright a-team source in LOFAR HBA
# they need to be in "./mss/"

import sys, os, glob, re
from shutil import copy2, copytree, move
import numpy as np
import casacore.tables as pt
import lsmtool

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-virgo.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-virgo.walker')

parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_virgo','parset_dir')
init_model = parset.get('LOFAR_virgo','init_model')
userReg = parset.get('model','userReg')
bl2flag = parset.get('flag','stations')

uvlambdamin = 30
t_int = 4
# nchan = 1
#############################################################################

field_model = '/beegfs/p1uy068/virgo/models/m87_field/'
# TODO
# Flag below 30 deg? Or can we use this?
# Amount of smoothing in BLSmooth?
# Solution intervalls? For FR probably we can go to at least 8 min?
# Station/antennaconstrains? For FR just all CS or maybe only superterp or so? For phases any constraint?
# updateweights
# amp normalization -> fulljones???

# Clear
with w.if_todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('img')
    os.makedirs('img')
    # here images, models, solutions for each group will be saved
    lib_util.check_rm('self')
    if not os.path.exists('self/solutions'): os.makedirs('self/solutions')
    if not os.path.exists('self/plots'): os.makedirs('self/plots')
    if not os.path.exists('self/masks'): os.makedirs('self/masks')
### DONE

MSs = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s )
if not MSs.isHBA:
    logger.error('Only HBA measurement sets supported.')
    sys.exit(1)
try:
    MSs.print_HAcov()
except:
    logger.error('Problem with HAcov, continue anyway.')

with pt.table(MSs.getListObj()[0].pathMS + '/OBSERVATION', ack=False) as t:
    if 'M87' in t.getcell('LOFAR_TARGET',0):
        target = 'VirA'
    else:
        raise ValueError('Target is not M87. Only M87 is implemented.')

basemask = 'self/masks/basemask.fits'
basemaskC = 'self/masks/basemaskC.fits'
baseregion = 'self/masks/baseregion.reg'
baseregionC = 'self/masks/baseregionC.reg'
phasecentre = MSs.getListObj()[0].getPhaseCentre()

#################################################################################################
# Initial flagging
with w.if_todo('flag'):
    logger.info('Flagging...')
    MSs.run(f'DP3 {parset_dir}/DP3-flag.parset msin=$pathMS ant.baseline=\"{bl2flag}\" '
            f'aoflagger.strategy={parset_dir}/HBAdefaultwideband.lua uvmin.uvlambdamin={uvlambdamin}',
            log='$nameMS_flag.log', commandType='DP3')

# Inital average -> for now to 8, 1 chan -> might wanna increase res later!
# Also cut frequencies above 168 MHz such that the number of channels is multiple of 4
if not os.path.exists('mss-avg'):
    timeint_init = MSs.getListObj()[0].getTimeInt()
    avgtimeint = int(round(t_int/timeint_init))  # to 16 seconds
    nchan_init = len(MSs.getFreqs())
    nchan = np.sum(np.array(MSs.getFreqs()) < 168e6) # only use 120-168 MHz
    nchan = 96 * (nchan // 96) # multiple of 48 after average
    os.makedirs('mss-avg')
    logger.info('Averaging in time (%is -> %is), channels: %ich -> %ich)' % (timeint_init,timeint_init*avgtimeint,nchan_init,nchan/2))
    MSs.run('DP3 '+parset_dir+'/DP3-avg.parset msin=$pathMS msout=mss-avg/$nameMS.MS msin.datacolumn=DATA msin.nchan='+str(nchan)+' \
            avg.timestep='+str(avgtimeint)+' avg.freqstep=2',
            log='$nameMS_initavg.log', commandType='DP3')
    # logger.info('Backup flags...')
    # MSs.addcol('FLAG_BKP', 'FLAG')
    # MSs.addcol('FLAG_ROW_BKP', 'FLAG_ROW')
    # MSs.run('taql "UPDATE $pathMS SET FLAG_BKP=FLAG"', log='$nameMS_taql.log', commandType='general')
    # MSs.run('taql "UPDATE $pathMS SET FLAG_ROW_BKP=FLAG_BKP"', log='$nameMS_taql.log', commandType='general')

MSs = lib_ms.AllMSs( glob.glob('mss-avg/TC*[0-9].MS'), s )

# Add model to MODEL_DATA
with w.if_todo('init_model'):
    if os.path.exists('model/img-0000-model.fits'):
        n = len(glob.glob('model/img-[0-9]*-model.fits'))
        logger.info('Predict (wsclean: %s - chan: %i)...' % ('model', n))
        s.add('wsclean -predict -name model/img -j ' + str(s.max_processors) + ' -channels-out ' + str(
            n) + ' ' + MSs.getStrWsclean(), \
              log='wscleanPRE-init.log', commandType='wsclean', processors='max')
        s.run(check=True)
    else:
        # copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
        sourcedb_basename = init_model.split('/')[-1]
        for MS in MSs.getListStr():
            lib_util.check_rm(MS + '/' + sourcedb_basename)
            logger.debug('Copy: ' + init_model + ' -> ' + MS)
            copy2(init_model, MS)

        # note: do not add MODEL_DATA or the beam is transported from DATA, while we want it without beam applied
        logger.info('Predict (DP3: %s - %s))...' % (sourcedb_basename, target))
        MSs.run(f'DP3 {parset_dir}/DP3-predict.parset msin=$pathMS pre.usebeammodel=false pre.sources={target} '
                f'pre.sourcedb=$pathMS/{sourcedb_basename}', log='$nameMS_pre.log', commandType='DP3')

#####################################################################################################
# Initial solve: FR
with w.if_todo('solve_fr'):
    logger.info('Add column CIRC_PHASEDIFF_DATA...')
    MSs.addcol('CIRC_PHASEDIFF_DATA', 'DATA', usedysco=False)
    logger.info('BL-smooth...')
    MSs.run('BLsmooth.py -r -c 4 -n 8 -f 5e-3 -i CIRC_PHASEDIFF_DATA -o CIRC_PHASEDIFF_DATA $pathMS',
            log='$nameMS_smooth.log', commandType='python', maxThreads=8)

    logger.info('Converting to circular...')
    MSs.run('mslin2circ.py -s -i $pathMS:CIRC_PHASEDIFF_DATA -o $pathMS:CIRC_PHASEDIFF_DATA',
            log='$nameMS_lincirc.log', commandType='python', maxThreads=2)

    # Get circular phase diff CIRC_PHASEDIFF_DATA -> CIRC_PHASEDIFF_DATA
    logger.info('Get circular phase difference...')
    MSs.run('taql "UPDATE $pathMS SET\
     CIRC_PHASEDIFF_DATA[,0]=0.5*EXP(1.0i*(PHASE(CIRC_PHASEDIFF_DATA[,0])-PHASE(CIRC_PHASEDIFF_DATA[,3]))), \
     CIRC_PHASEDIFF_DATA[,3]=CIRC_PHASEDIFF_DATA[,0], \
     CIRC_PHASEDIFF_DATA[,1]=0+0i, \
     CIRC_PHASEDIFF_DATA[,2]=0+0i"', log='$nameMS_taql_phdiff.log', commandType='general')

    logger.info('Creating FR_MODEL_DATA...')  # take from MODEL_DATA but overwrite
    MSs.addcol('FR_MODEL_DATA', 'MODEL_DATA', usedysco=False)
    MSs.run('taql "UPDATE $pathMS SET FR_MODEL_DATA[,0]=0.5+0i, FR_MODEL_DATA[,1]=0.0+0i, FR_MODEL_DATA[,2]=0.0+0i, \
     FR_MODEL_DATA[,3]=0.5+0i"', log='$nameMS_taql_frmodel.log', commandType='general')

    # Solve cal_SB.MS:CIRC_PHASEDIFF_DATA against FR_MODEL_DATA (only solve)
    logger.info('Solving circ phase difference ...')
    MSs.run(f'DP3 {parset_dir}/DP3-solFR.parset msin=$pathMS sol.h5parm=$pathMS/fr.h5 sol.uvlambdamin={uvlambdamin} sol.solint={120//t_int}',
            log='$nameMS_solFR.log', commandType="DP3")
    lib_util.run_losoto(s, f'fr', [ms + '/fr.h5' for ms in MSs.getListStr()],
                        [parset_dir + '/losoto-fr.parset'])
    move('cal-fr.h5', 'self/solutions/')
    move('plots-fr', 'self/plots/')
    # Delete cols again to not waste space
    MSs.run('taql "ALTER TABLE $pathMS DELETE COLUMN CIRC_PHASEDIFF_DATA, FR_MODEL_DATA"',
            log='$nameMS_taql_delcol.log', commandType='general')
### DONE

with w.if_todo('cor_fr'):
    # Correct FR  DATA -> FR_CORRECTED_DATA
    # MSs.addcol('CIRC_PHASEDIFF_DATA', 'CORRECTED_DATA', usedysco=False)
    logger.info('Correcting FR...')
    MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=FR_CORRECTED_DATA \
            cor.parmdb=self/solutions/cal-fr.h5 cor.correction=rotationmeasure000', log='$nameMS_corFR.log',
            commandType='DP3')
    # BL-smooth FR_CORRECTED_DATA -> FR_SMOOTHED_DATA
    logger.info('BL-smooth...')
    MSs.run('BLsmooth.py -r -c 4 -n 8 -f 5e-3 -i FR_CORRECTED_DATA -o FR_SMOOTHED_DATA $pathMS',
            log='$nameMS_smooth.log', commandType='python', maxThreads=8)
### DONE
# Self-cal cycle
field_subtracted = False
for c in range(100):
    # Solve cal_SB.MS: FR_SMOOTHED_DATA (only solve)
    with w.if_todo('solve_iono_c%02i' % c):
        logger.info('Solving scalarphase...')
        # TODO sol.smoothnessconstraint=1e6 FR_SMOOTHED_DATA
        MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=FR_CORRECTED_DATA '
                f'sol.h5parm=$pathMS/iono.h5 sol.mode=scalarcomplexgain  sol.uvlambdamin={uvlambdamin}' , log=f'$nameMS_sol_iono-c{c}.log', commandType="DP3")

        lib_util.run_losoto(s, f'iono-c{c:02}', [ms+'/iono.h5' for ms in MSs.getListStr()], \
                            [#parset_dir+'/losoto-flag.parset',
                             parset_dir+'/losoto-plot-scalaramp.parset',
                             parset_dir+'/losoto-plot-scalarph.parset'])

        move(f'cal-iono-c{c:02}.h5', 'self/solutions/')
        move(f'plots-iono-c{c:02}', 'self/plots/')
        # Correct all FR_CORRECTED_DATA -> CORRECTED_DATA

    with w.if_todo('cor_iono_c%02i' % c):
        logger.info('Scalarphase correction...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=FR_CORRECTED_DATA cor.updateweights=False cor.parmdb=self/solutions/cal-iono-c{c:02}.h5 cor.correction=phase000', \
                log=f'$nameMS_cor_iono-c{c:02}.log', commandType='DP3')
        ### DONE


    if 0 < c < 4 :
        # Solve cal_SB.MS: CORRECTED_DATA --smooth--> SMOOTHED_DATA --solve-->
        with w.if_todo('solve_gain_c%02i' % c):
            logger.info('BL-smooth...')
            MSs.run('BLsmooth.py -r -c 4 -n 8 -f 5e-3 -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS',
                    log='$nameMS_smooth.log', commandType='python', maxThreads=8)
            logger.info('Solving gain...')
            MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA '
                    f'sol.h5parm=$pathMS/gain.h5 sol.mode=diagonal sol.smoothnessconstraint=2e6 sol.nchan=6 sol.solint={120//t_int} '
                    f'sol.uvlambdamin={uvlambdamin}', log=f'$nameMS_sol_gain-c{c}.log',
                    commandType="DP3")

            lib_util.run_losoto(s, f'gain-c{c:02}', [ms + '/gain.h5' for ms in MSs.getListStr()], \
                                [parset_dir + '/losoto-flag.parset',
                                 parset_dir + '/losoto-norm.parset',
                                 parset_dir + '/losoto-plot-amp.parset',
                                 parset_dir + '/losoto-plot-ph.parset'])

            move(f'cal-gain-c{c:02}.h5', 'self/solutions/')
            move(f'plots-gain-c{c:02}', 'self/plots/')

        # Correct gain amp and ph CORRECTED_DATA -> CORRECTED_DATA
        with w.if_todo('cor_gain_c%02i' % c):
            logger.info('Gain correction...')
            MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.steps=[amp,phase] \
                    cor.amp.parmdb=self/solutions/cal-gain-c{c:02}.h5 cor.amp.correction=amplitude000 '
                    f'cor.parmdb=self/solutions/cal-gain-c{c:02}.h5 cor.phase.correction=phase000',
                    log=f'$nameMS_cor_gain-c{c:02}.log',
                    commandType='DP3')
            ### DONE
    elif c >= 4:
        with w.if_todo('solve_fulljones%02i' % c):
            logger.info('BL-smooth...')
            MSs.run('BLsmooth.py -r -c 4 -n 8 -f 5e-3 -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS',
                    log='$nameMS_smooth.log', commandType='python', maxThreads=8)
            logger.info('Solving full-Jones...')
            # solint from 180 to 60, chan from 4 to 1, smoothness from 3 to 2
            MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA '
                    f'sol.h5parm=$pathMS/fulljones.h5 sol.mode=fulljones sol.smoothnessconstraint=2e6 sol.nchan=1 sol.solint={60//t_int} '
                    f'sol.uvlambdamin={uvlambdamin}', log=f'$nameMS_sol_fulljones-c{c}.log',
                    commandType="DP3")

            lib_util.run_losoto(s, f'fulljones-c{c:02}', [ms + '/fulljones.h5' for ms in MSs.getListStr()], \
                                [#parset_dir + '/losoto-norm.parset',
                                 parset_dir + '/losoto-plot-fulljones.parset'])

            move(f'cal-fulljones-c{c:02}.h5', 'self/solutions/')
            move(f'plots-fulljones-c{c:02}', 'self/plots/')

        # Correct gain amp and ph CORRECTED_DATA -> CORRECTED_DATA
        with w.if_todo('cor_fulljones_c%02i' % c):
            logger.info('Full-Jones correction...')
            MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA '
                    f'cor.correction=fulljones cor.parmdb=self/solutions/cal-fulljones-c{c:02}.h5 '
                    f'cor.soltab=\[amplitude000,phase000\]', log=f'$nameMS_cor_gain-c{c:02}.log', commandType='DP3')

    ###################################################################################################################
    # clean CORRECTED_DATA
    imagename = f'img/img-c{c:02}'
    wsclean_params = {
        'scale': '0.5arcsec',
        'size': 2400,
        'weight': 'briggs -1.5',
        'join_channels': '',
        # 'deconvolution_channels': 32,
        'fit_spectral_pol': 8,  # 3 worked fine, let's see if the central residual improves with 5
        'channels_out': len(MSs.getFreqs()) // 12,
        'minuv_l': uvlambdamin,
        'multiscale': '',
        'name': imagename,
        'no_update_model_required': '',
        'do_predict': True,
        'baseline_averaging': 10
    }
    with w.if_todo('imaging_c%02i' % c):
        logger.info('Cleaning (cycle: '+str(c)+')...')

        if not os.path.exists(basemask) or not os.path.exists(basemaskC):
            logger.info('Create masks...')
            # dummy clean to get image -> mask
            lib_util.run_wsclean(s, 'wsclean-c' + str(c) + '.log', MSs.getStrWsclean(), niter=0, channel_range='0 1',
                                 interval='0 10', name=imagename, scale=wsclean_params['scale'], size=wsclean_params['size'], nmiter=0)
            # create basemask
            copy2(f'{parset_dir}/masks/VirAhbaEllipse.reg', f'{baseregion}')
            copy2(f'{imagename}-image.fits', f'{basemask}')
            lib_img.blank_image_reg(basemask, baseregion, inverse=True, blankval=0.)
            lib_img.blank_image_reg(basemask, baseregion, inverse=False, blankval=1.)
            if userReg != '':
                lib_img.blank_image_reg(basemask, userReg, inverse=True, blankval=1.)

        logger.info('Cleaning...')
        lib_util.run_wsclean(s, f'wsclean-c{c}.log', MSs.getStrWsclean(), niter=1500000,
                             fits_mask=basemask, multiscale_scales='0,20,30,45,66,99,150', nmiter=30, mgain=0.6, gain=0.08, multiscale_gain=0.12,
                             auto_threshold=1.2, auto_mask=3.0, **wsclean_params)
        os.system(f'cat logs/wsclean-c{c}.log | grep "background noise"')

    # widefield_model = False
    # if c >= 5 and not field_subtracted:
    #     with w.if_todo('imaging_wide_c%02i' % c):
    #         logger.info('SET SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA')
    #         MSs.addcol('SUBTRACTED_DATA', 'CORRECTED_DATA')
    #         MSs.run('taql "UPDATE $pathMS SET SUBTRACTED_DATA = CORRECTED_DATA-MODEL_DATA"', log='$nameMS_taql_subtract.log', commandType='general')
    #
    #         logger.info('Cleaning Virgo A subtracted wide-field image...')
    #         lib_util.run_wsclean(s, f'wsclean-wide-c{c}.log', MSs.getStrWsclean(), weight='briggs -0.5', data_column='SUBTRACTED_DATA',
    #                              name=imagename+'-wide', parallel_deconvolution=1024, scale='2.0arcsec', size=7200, niter=500000,
    #                              join_channels='', channels_out=6, nmiter=15, fit_spectral_pol=3, minuv_l=uvlambdamin, multiscale='', multiscale_max_scales=5,
    #                              mgain=0.85, auto_threshold=1.0, auto_mask=4.0, baseline_averaging='', no_update_model_required='', do_predict=True, local_rms='',
    #                              fits_mask=parset_dir+'/masks/FieldMask.fits')
    #         widefield_model = True
    #         os.system(f'cat logs/wsclean-wide-c{c}.log | grep "background noise"')

            # logger.info('Cleaning Virgo A subtracted low-res wide-field image...')
            # lib_util.run_wsclean(s, f'wsclean-wide-lr-c{c}.log', MSs.getStrWsclean(), weight='briggs -0.5', taper_gaussian='30arcsec', data_column='SUBTRACTED_DATA',
            #                      name=imagename+'-wide-lr', parallel_deconvolution=2048, scale='6.0arcsec', size=3000, niter=500000,
            #                      join_channels='', channels_out=6, fit_spectral_pol=3, minuv_l=uvlambdamin, multiscale='', multiscale_max_scales=5,
            #                      mgain=0.85, auto_threshold=1.0, auto_mask=3.0, baseline_averaging='', no_update_model_required='')
            # os.system(f'cat logs/wsclean-wide-c{c}.log | grep "background noise"')

        # with w.if_todo('subtract_wide_c%02i' % c):
            # if not widefield_model:
            #     logger.info('Predict rest-field...')
            #     n = len(glob.glob(field_model + 'm87-field-[0-9]*-model.fits'))
            #     logger.info('Predict (wsclean: %s - chan: %i)...' % ('model-field', n))
            #     s.add(f'wsclean -predict -name {field_model}m87-field -j {s.max_processors} -channels-out {n} {MSs.getStrWsclean()}',
            #           log='wscleanPRE-field.log', commandType='wsclean', processors='max')
            #     s.run(check=True)
            #
            # logger.info('TEST EMPTY')
            # logger.info('Set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA')
            # MSs.run('taql "UPDATE $pathMS SET SUBTRACTED_DATA = SUBTRACTED_DATA-MODEL_DATA"', log='$nameMS_taql_subtract_empty.log',
            #         commandType='general')
            # lib_util.run_wsclean(s, f'wsclean-empty-c{c}.log', MSs.getStrWsclean(), weight='briggs -0.5', data_column='SUBTRACTED_DATA',
            #                      name=imagename+'-empty', scale='2.0arcsec', size=7200, minuv_l=uvlambdamin, no_update_model_required='')

            # logger.info('Corrupt widefield MODEL_DATA...')
            # logger.info('Scalarphase corruption (MODEL_DATA)...')
            # MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA '
            #         f'cor.updateweights=False cor.parmdb=self/solutions/cal-iono-c{c:02}.h5 cor.correction=phase000 '
            #         f'cor.invert=False', log=f'$nameMS_corrupt_iono-c{c:02}.log', commandType='DP3')
            # logger.info('Full-Jones corruption (MODEL_DATA)...')
            # MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA '
            #         f'cor.correction=fulljones cor.parmdb=self/solutions/cal-fulljones-c{c:02}.h5 '
            #         f'cor.soltab=\[amplitude000,phase000\] cor.invert=False', log=f'$nameMS_corrupt_gain-c{c:02}.log', commandType='DP3')
            #
            # logger.info('Set FR_CORRECTED_DATA = FR_CORRECTED_DATA - MODEL_DATA')
            # MSs.run('taql "UPDATE $pathMS SET FR_CORRECTED_DATA = FR_CORRECTED_DATA-MODEL_DATA"', log='$nameMS_taql_subtract.log',
            #         commandType='general')
            #
            # logger.info('BL-smooth...')
            # MSs.run('BLsmooth.py -r -c 4 -n 8 -f 5e-3 -i FR_CORRECTED_DATA -o FR_SMOOTHED_DATA $pathMS',
            #         log='$nameMS_smooth.log', commandType='python', maxThreads=8)
            #
            # logger.info('Get back Virgo A MODEL_DATA...')
            # s.add(f'wsclean -predict -name {imagename} -j {s.max_processors} -channels-out {wsclean_params["channels_out"]} {MSs.getStrWsclean()}',
            #       log='wscleanPRE-field.log', commandType='wsclean', processors='max')
            # s.run(check=True)
logger.info("Done.")
