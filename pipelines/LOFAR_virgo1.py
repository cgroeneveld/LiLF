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
import astropy.io.fits as pyfits

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-virgo1.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-virgo1.walker')

parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_virgo','parset_dir')
init_model = parset.get('model','fits_model')
userReg = parset.get('model','userReg')
bl2flag = parset.get('flag','stations')

#############################################################################

# m87_model = '/beegfs/p1uy068/virgo/models/m87-9/m87'
# TODO
# Solution intervals? For FR probably we can go to at least 8 min?
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
    os.makedirs('self/solutions')
    os.makedirs('self/plots')
    os.makedirs('self/masks')
### DONE

uvlambdamin = 30
uvwmmax = 90e3 # this will remove RS210-RS508 and RS210-RS509
t_int = 4

basemask = 'self/masks/basemask.fits'
basemaskC = 'self/masks/basemaskC.fits'
baseregion = 'self/masks/baseregion.reg'
baseregionC = 'self/masks/baseregionC.reg'

#################################################################################################
# Initial flagging
MSs = lib_ms.AllMSs( glob.glob('mss-slim/TC*[0-9].MS'), s )
phasecentre = MSs.getListObj()[0].getPhaseCentre()
#
# Add model to MODEL_DATA
# with w.if_todo('init_model'):
#     logger.info('Predict (wsclean: %s )...' % (m87_model))
#     s.add(f'wsclean -predict -name {m87_model} -j {s.max_processors} {MSs.getStrWsclean()}',
#           log='wscleanPRE-init.log', commandType='wsclean', processors='max')
#     s.run(check=True)

#####################################################################################################
### DONE


# Self-cal cycle
field_subtracted = False
for c in range(100):
    # Solve cal_SB.MS: FR_SMOOTHED_DATA (only solve)
    with w.if_todo('solve_iono_c%02i' % c):
        logger.info('Solving scalarphase...')
        MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=DATA '
                f'sol.h5parm=$pathMS/iono.h5 sol.mode=scalarcomplexgain sol.smoothnessconstraint=1e6 '
                f'sol.uvlambdamin={uvlambdamin} sol.uvmmax={uvwmmax}' , log=f'$nameMS_sol_iono-c{c}.log',
                commandType="DP3")

        lib_util.run_losoto(s, f'iono-c{c:02}', [ms+'/iono.h5' for ms in MSs.getListStr()], \
                            [#parset_dir+'/losoto-flag.parset',
                             parset_dir+'/losoto-plot-scalaramp.parset',
                             parset_dir+'/losoto-plot-scalarph.parset'])

        move(f'cal-iono-c{c:02}.h5', 'self/solutions/')
        move(f'plots-iono-c{c:02}', 'self/plots/')
        # Correct all FR_CORRECTED_DATA -> CORRECTED_DATA

    with w.if_todo('cor_iono_c%02i' % c):
        logger.info('Scalarphase correction...')
        MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA cor.updateweights=False cor.parmdb=self/solutions/cal-iono-c{c:02}.h5 cor.correction=phase000', \
                log=f'$nameMS_cor_iono-c{c:02}.log', commandType='DP3')
        ### DONE

    if False:
        # Solve cal_SB.MS: CORRECTED_DATA --smooth--> SMOOTHED_DATA --solve-->
        with w.if_todo('solve_gain_c%02i' % c):
            logger.info('Solving gain...')
            MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA '
                    f'sol.h5parm=$pathMS/gain.h5 sol.mode=diagonal sol.smoothnessconstraint=4e6 sol.nchan=4 sol.solint={900//t_int} '
                    f'sol.uvlambdamin={uvlambdamin} sol.uvmmax={uvwmmax}', log=f'$nameMS_sol_gain-c{c}.log',
                    commandType="DP3")

            lib_util.run_losoto(s, f'gain-c{c:02}', [ms + '/gain.h5' for ms in MSs.getListStr()], \
                                [parset_dir + '/losoto-flag.parset',
                                 parset_dir + '/losoto-diag.parset'])

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
    else:
        with w.if_todo('solve_gain_c%02i' % c):
            logger.info('Solving full-Jones...')
            # solchan = len(MSs.getFreqs())//(len(MSs.getFreqs()) // 12)
            MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA '
                    f'sol.h5parm=$pathMS/fulljones.h5 sol.mode=fulljones sol.nchan=4 sol.solint={300//t_int} '
                    f'sol.uvlambdamin={uvlambdamin} sol.smoothnessconstraint=4e6', log=f'$nameMS_sol_fulljones-c{c}.log',
                    commandType="DP3")

            lib_util.run_losoto(s, f'fulljones-c{c:02}', [ms + '/fulljones.h5' for ms in MSs.getListStr()], \
                                [#parset_dir + '/losoto-norm.parset',
                                 parset_dir + '/losoto-plot-fulljones.parset'])

            move(f'cal-fulljones-c{c:02}.h5', 'self/solutions/')
            move(f'plots-fulljones-c{c:02}', 'self/plots/')

        # Correct gain amp and ph CORRECTED_DATA -> CORRECTED_DATA
        with w.if_todo('cor_gain_c%02i' % c):
                logger.info('Full-Jones correction...')
                MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA '
                        f'cor.correction=fulljones cor.parmdb=self/solutions/cal-fulljones-c{c:02}.h5 '
                        f'cor.soltab=\[amplitude000,phase000\]', log=f'$nameMS_cor_gain-c{c:02}.log', commandType='DP3')
    #
    ###################################################################################################################
    # clean CORRECTED_DATA
    imagename = f'img/img-c{c:02}'
    # TODO test NO fit spectrum + NO join channels
    wsclean_params = {
        'scale': f'{MSs.resolution/5}arcsec', #'0.5arcsec',
        'size': int(1600/(MSs.resolution/5)),
        'weight': 'briggs -0.8', # IS?
        'minuv_l': uvlambdamin,
        'name': imagename,
        'no_update_model_required': '',
        'do_predict': True,
        'mgain': 0.85,
    }
    with w.if_todo('imaging_c%02i' % c):
        logger.info(f'Cleaning (cycle: {c}; size: {wsclean_params["size"]}pix scale: {wsclean_params["scale"]})...')

        if not os.path.exists(basemask) or not os.path.exists(basemaskC):
            logger.info('Create masks...')
            # dummy clean to get image -> mask
            lib_util.run_wsclean(s, 'wsclean-c' + str(c) + '.log', MSs.getStrWsclean(), niter=0, channel_range='0 1',
                                 interval='0 10', name=imagename, scale=wsclean_params['scale'], size=wsclean_params['size'], nmiter=0)
            # create basemask
            copy2(f'{parset_dir}/masks/VirAhba.reg', f'{baseregion}')
            copy2(f'{imagename}-image.fits', f'{basemask}')
            copy2(f'{parset_dir}/masks/VirAChba.reg', f'{baseregionC}')
            copy2(f'{imagename}-image.fits', f'{basemaskC}')
            lib_img.blank_image_reg(basemask, baseregion, inverse=True, blankval=0.)
            lib_img.blank_image_reg(basemask, baseregion, inverse=False, blankval=1.)
            lib_img.blank_image_reg(basemaskC, baseregionC, inverse=True, blankval=0.)
            lib_img.blank_image_reg(basemaskC, baseregionC, inverse=False, blankval=1.)
            if userReg != '':
                lib_img.blank_image_reg(basemask, userReg, inverse=True, blankval=1.)
        # iterate cocoon - halo
        switch_gain = 0.77
        threshold = 0.1 # initial threshold
        threshold_final = 0.002
        last_it = False
        i = 0
        am = dict()
        while threshold >= threshold_final:
            logger.info(f'Iteration threshold = {threshold}')
            # Clean delta scale on cocoon
            logger.info('Cleaning cocoon...')
            lib_util.run_wsclean(s, f'wsclean1-c{c}.log', MSs.getStrWsclean(), niter=1000000,
                                 fits_mask=basemaskC, gain=0.05, threshold=threshold, **wsclean_params)
            os.system(f'cat logs/wsclean1-c{c}.log | grep "background noise"')
            try:
                del wsclean_params['mgain']  # only use C-S during first call where many orders of CC are found
            except KeyError:
                pass
            wsclean_params['cont'] = True
            wsclean_params['reuse-psf'] = imagename
            with pyfits.open(imagename + '-residual.fits') as fits: # MFS
                peak_halo_resid = np.max(np.abs(fits[0].data))
                print(f'Peak residual in halo: {peak_halo_resid}')

            if threshold < 4 * threshold_final:
                am['auto_mask'] = 3.0
            if peak_halo_resid > threshold:
                logger.info('Cleaning (multi-scale) halo...')
                lib_util.run_wsclean(s, f'wsclean2-c{c}.log', MSs.getStrWsclean(), niter=100000, multiscale='', multiscale_scales='0,20,40,80,160,320',
                                     fits_mask=basemask, threshold=threshold,  **{**wsclean_params, **am})
                os.system(f'cat logs/wsclean2-c{c}.log | grep "background noise"')
            i += 1
            threshold = (1.0 - switch_gain) * threshold
            if last_it:
                break
            elif threshold < threshold_final: # perform last iteration using the minimal threshold
                threshold = threshold_final
                last_it = True

logger.info("Done.")
