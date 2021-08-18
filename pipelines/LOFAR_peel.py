#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Peel a bright source outside of the PB from the data
# Two possible ways to run this pipeline:
# 1. You provide either a ds9 circle/polygon peelReg (default 'peel.reg') and either a sourcedb or a fits_model
#    in the [model] section. The pipeline will solve against this initial model and then continue w/ self-cal.
# The MS files need to be in "./mss/"

import sys, os, glob, re
from shutil import copy2, copytree, move
import numpy as np
import astropy.io.fits as fits
import casacore.tables as pt
import lsmtool

from LiLF import lib_ms, lib_img, lib_util, lib_log
#############################################################################
logger_obj = lib_log.Logger('pipeline-peel.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-peel.walker')

parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_peel','parset_dir')
cal_dir = parset.get('LOFAR_peel','cal_dir')
data_dir = parset.get('LOFAR_peel','data_dir')
peelReg = parset.get('LOFAR_peel','peelReg')  # default peel.reg
predictReg = parset.get('LOFAR_peel','predictReg')  # default None
sourcedb = parset.get('model','sourcedb')
fits_model = parset.get('model','fits_model')
bl2flag = parset.get('flag','stations')

# Some script-wide definitions
uvlambdamin = 30 # why not 100? What is used in ddf-pipeline?
t_int = 8
final_chan_per_sb = 2
pixscale = '1.5arcsec'
imgsize = 1600
#############################################################################

def solve_and_apply(MSs_object, suffix, sol_factor_t=1, sol_factor_f=1):
    """
    Run a scalarphase-solve and a fulljones-solve on MSs_object and apply the solutions to the data each time.
    Parameters
    ----------
    MSs_object
    suffix: string, suffix for the logs, i.e. 'c00' will give us xxx-c00.log
    sol_factor_t
    sol_factor_f
    """
    with w.if_todo(f'solve_iono_{suffix}'):
        logger.info('Solving scalarphase...')
        s.add(f'DP3 {parset_dir}/DP3-soldd.parset msin=[{",".join(MSs_object.getListStr())}] sol.h5parm=$pathMS/iono.h5 '
                         f'sol.mode=scalarcomplexgain sol.nchan={sol_factor_f:d} sol.smoothnessconstraint={1.5e6} '
                         f'sol.solint={sol_factor_t:d} sol.uvlambdamin={uvlambdamin}',
                         log=f'sol_iono-{suffix}.log', commandType="DP3")
        s.run(check=True)

        lib_util.run_losoto(s, f'iono-c{c:02}', [ms + '/iono.h5' for ms in MSs_object.getListStr()], \
                            [  # parset_dir+'/losoto-flag.parset',
                                parset_dir + '/losoto-plot-scalaramp.parset',
                                parset_dir + '/losoto-plot-scalarph.parset'])
        move(f'cal-iono-{suffix}.h5', 'peel/solutions/')
        move(f'plots-iono-{suffix}', 'peel/plots/')

    with w.if_todo(f'cor_iono_{suffix}'):
        logger.info('Scalarphase correction DATA -> CORRECTED_DATA...')
        MSs_object.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA cor.updateweights=False '
                         f'cor.parmdb=peel/solutions/cal-iono-{suffix}.h5 cor.correction=phase000', \
                         log=f'$nameMS_cor_iono-{suffix}.log', commandType='DP3')
        ### DONE

    with w.if_todo(f'solve_fulljones_{suffix}'):
        logger.info('Solving full-Jones...')
        s.add(f'DP3 {parset_dir}/DP3-soldd.parset msin=[{",".join(MSs_object.getListStr())}] msin.datacolumn=CORRECTED_DATA '
                         f'sol.h5parm=$pathMS/fulljones.h5 sol.mode=fulljones sol.nchan={4 * sol_factor_f:d} sol.solint={sol_factor_t * 128 // t_int:d} '
                         f'sol.smoothnessconstraint={2.0e6} sol.uvlambdamin={uvlambdamin}',
                         log=f'sol_fulljones-{suffix}.log', commandType="DP3")
        s.run(check=True)

        lib_util.run_losoto(s, f'fulljones-{suffix}', [ms + '/fulljones.h5' for ms in MSs_object.getListStr()], \
                            [parset_dir + '/losoto-plot-fulljones.parset'])
        move(f'cal-fulljones-{suffix}.h5', 'peel/solutions/')
        move(f'plots-fulljones-{suffix}', 'peel/plots/')

    # Correct gain amp and ph CORRECTED_DATA -> CORRECTED_DATA
    with w.if_todo(f'cor_fulljones_{suffix}'):
        logger.info('Full-Jones correction CORRECTED_DATA -> CORRECTED_DATA...')
        MSs_object.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.correction=fulljones '
                         f'cor.parmdb=peel/solutions/cal-fulljones-{suffix}.h5 cor.soltab=\[amplitude000,phase000\]',
                         log=f'$nameMS_cor_gain-{suffix}.log', commandType='DP3')

def corrupt_subtract_testimage(MSs_object, suffix, sol_suffix=None):
    """
    Use solutions of a suffix to corrupt and subtract MODEL_DATA, also create an empty test image.
    Parameters
    ----------
    MSs_object
    suffix
    """
    if sol_suffix is None:
        sol_suffix = suffix
    with w.if_todo(f'subtract-{suffix}'):
        logger.info('Scalarphase corruption...')
        MSs_object.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA '
                      f'msout.datacolumn=CORRUPTED_MODEL_DATA cor.updateweights=False '
                      f'cor.parmdb=peel/solutions/cal-iono-{sol_suffix}.h5 cor.correction=phase000 cor.invert=False',
                      log=f'$nameMS_corrup_iono-{sol_suffix}.log', commandType='DP3')

        logger.info('Full-Jones corruption...')
        MSs_object.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRUPTED_MODEL_DATA '
                      f'msout.datacolumn=CORRUPTED_MODEL_DATA cor.correction=fulljones '
                      f'cor.parmdb=peel/solutions/cal-fulljones-{sol_suffix}.h5 cor.soltab=\[amplitude000,phase000\] '
                      f'cor.invert=False', log=f'$nameMS_corrup_gain-{sol_suffix}.log', commandType='DP3')

        logger.info('SET SUBTRACTED_DATA = DATA - CORRUPTED_MODEL_DATA')
        MSs_object.addcol('SUBTRACTED_DATA', 'DATA')
        MSs_object.run('taql "UPDATE $pathMS SET SUBTRACTED_DATA = DATA - CORRUPTED_MODEL_DATA"',
                      log='$nameMS_taql_subtract.log',
                      commandType='general')

    with w.if_todo(f'correct-subtracted-{suffix}'):
        logger.info(f'Scalarphase correction...')
        MSs_object.run(
            f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBTRACTED_DATA '
            f'msout.datacolumn=SUBTRACTED_DATA cor.updateweights=False '
            f'cor.parmdb=peel/solutions/cal-iono-{sol_suffix}.h5 cor.correction=phase000', \
            log=f'$nameMS_cor_iono-{sol_suffix}.log', commandType='DP3')

    with w.if_todo(f'test-image-empty-{suffix}'):
        logger.info('Test empty... (SUBTRACTED_DATA)')
        lib_util.run_wsclean(s, f'wsclean-peel.log', MSs_object.getStrWsclean(), weight='briggs -0.5',
                             data_column='SUBTRACTED_DATA',
                             name=imagename + '-empty', scale='2.0arcsec', size=2000, niter=0, nmiter=0,
                             no_update_model_required='',
                             minuv_l=uvlambdamin)

def predict_fits_model(MSs_object, model_basename, stepname='init_model', predict_reg=None):
    """
    Predict a fits_model into one or more MS files. Optionally, blank fits model with 'predict_reg'
    Parameters
    ----------
    MSs_object
    model_basename
    stepname
    predict_reg
    """
    with w.if_todo(stepname):
        if predict_reg:
            for model_img in glob.glob(model_basename + '*-model.fits'):
                lib_img.blank_image_reg(model_img, predict_reg, inverse=True, blankval=0.)
        assert os.path.isfile(model_basename + '-0000-model.fits')
        logger.info(f'Using fits model {model_basename}.')
        n = len(glob.glob(model_basename + '-[0-9]*-model.fits'))
        logger.info('Predict (wsclean: %s - chan: %i)...' % ('model', n))
        s.add(f'wsclean -predict -name {model_basename}  -j {s.max_processors} -use-wgridder -channels-out {n} '
              f'{MSs_object.getStrWsclean()}', log=f'wscleanPRE-{stepname}.log', commandType='wsclean', processors='max')
        s.run(check=True)
#############################################################################
# Clear
with w.if_todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('img')
    os.makedirs('img')
    # here images, models, solutions for each group will be saved
    lib_util.check_rm('peel')
    os.makedirs('peel/solutions')
    os.makedirs('peel/plots')
    os.makedirs('peel/masks')
### DONE

MSs = lib_ms.AllMSs( glob.glob(data_dir + '/*.MS'), s, check_flags=False)
if not MSs.isHBA:
    logger.error('Only HBA measurement sets supported for now.')
    sys.exit(1)
try:
    MSs.print_HAcov()
except:
    logger.error('Problem with HAcov, continue anyway.')
logger.info(f'Resolution {MSs.resolution}\'\', using uvlambdamin={uvlambdamin}, averaging to {final_chan_per_sb}chan/SB in freq and to an '
            f'integration time of t_int={t_int}')

peelMask = 'peel/masks/peelmask.fits'
phasecentre = MSs.getListObj()[0].getPhaseCentre()
# region must be a list of ds9 circles and polygons (other shapes can be implemented in lib_util.Rgion_helper()
peelReg = lib_util.Region_helper(peelReg)
assert os.path.isfile(fits_model + '-0000-model.fits')
model_hdr = fits.open(glob.glob(fits_model + '-[0-9]*-model.fits')[0])[0].header
model_centre = np.array([model_hdr['CRVAL1'], model_hdr['CRVAL2']])
pointing_distance = lib_util.distanceOnSphere(*model_centre, *phasecentre)
logger.info(f"Distance between model and MS phase centre: {pointing_distance:.5f}deg")
##################################################
with w.if_todo('apply'):
    # Find solutions to apply
    cal_dirs = glob.glob(cal_dir + 'id*_-_3[c|C]196') + glob.glob(cal_dir + 'id*_-_3[c|C]295')
    if len(cal_dirs) == 0:
        logger.error(f'No calibrators found in cal dir: {cal_dir}')
        sys.exit(1)

    cal_times = []  # mean times of the cal
    for cal in cal_dirs:
        cal_ms = lib_ms.MS(glob.glob(f'{cal}/*.MS')[0]) # assume all MS are equal
        assert cal_ms.isCalibrator()
        cal_times.append(np.mean(cal_ms.getTimeRange()))
    obs_time = np.mean(MSs.getListObj()[0].getTimeRange())
    delta_t = np.abs(obs_time - np.array(cal_times))
    cal_id = np.argmin(delta_t)
    cal_dir =  cal_dirs[cal_id]
    if delta_t[cal_id] < 5*3600:
        logger.info(f'Found calibrator dir {cal_dirs[cal_id]} which is {delta_t[cal_id]/3600:.2f}h from mean observation time.')
    else:
        logger.error(f'Found calibrator dir {cal_dirs[cal_id]} which is {delta_t[cal_id]/3600:.2f}h from mean observation time!')

    logger.info('Calibrator directory: %s' % cal_dir)
    h5_pa = cal_dir + '/cal-pa.h5'
    h5_amp = cal_dir + '/cal-amp.h5'
    h5_iono = cal_dir + '/cal-iono.h5'
    if not os.path.exists(h5_pa) or not os.path.exists(h5_amp) or not os.path.exists(h5_iono):
        logger.error("Missing solutions in %s" % cal_dir)
        sys.exit()

    # Apply cal sol - SB.MS:DATA -> SB.MS:CORRECTED_DATA (polalign corrected)
    logger.info('Apply solutions (pa)...')
    MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA \
            cor.parmdb=' + h5_pa + ' cor.correction=polalign', log='$nameMS_cor1.log', commandType='DP3')

    # Apply cal sol - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (polalign corrected, calibrator corrected+reweight, beam corrected+reweight)
    logger.info('Apply solutions (amp/ph)...')
    MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.steps=[amp,clock] \
            cor.amp.parmdb=' + h5_amp + ' cor.amp.correction=amplitudeSmooth cor.amp.updateweights=True\
            cor.clock.parmdb=' + h5_iono + ' cor.clock.correction=clockMed000', log='$nameMS_cor2.log',
            commandType='DP3')

    # Beam correction CORRECTED_DATA -> CORRECTED_DATA (polalign corrected, beam corrected+reweight)
    logger.info('Beam correction...')
    MSs.run('DP3 ' + parset_dir + '/DP3-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_beam.log',
            commandType='DP3')
### DONE

#################################################################################################
# Initial flagging
with w.if_todo('flag'):
    logger.info('Flagging...')
    s.add(f'DP3 {parset_dir}/DP3-flag.parset msin=[{",".join(MSs.getListStr())}]  ant.baseline=\"{bl2flag}\" msin.datacolumn=CORRECTED_DATA '
            f'aoflagger.strategy={parset_dir}/HBAdefaultwideband.lua uvmin.uvlambdamin={uvlambdamin}',
            log='$nameMS_flag.log', commandType='DP3')
    s.run(check=True)
    logger.info('Remove bad timestamps...')
    MSs.run('flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')

    # need to make this work for more than one MS!
    # logger.info('Plot weights...')
    # MSs.run(f'reweight.py $pathMS -v -p -a "CS001HBA0"',
    #         log='$nameMS_weights.log', commandType='python')
    # os.system('move *.png self/plots')
### DONE
# Also cut frequencies above 168 MHz such that the number of channels is multiple of 4
with pt.table(MSs.getListObj()[0].pathMS + "/OBSERVATION") as tab:
    field = tab[0]["LOFAR_TARGET"][0]

if pointing_distance > 1/3600: # CASE 1 -> model and MS not aligned, peel
    if not os.path.exists('mss-shift'):
        timeint_init = MSs.getListObj()[0].getTimeInt()
        nchan_init = len(MSs.getFreqs())
        nchan = np.sum(np.array(MSs.getFreqs()) < 168.3e6)  # only use 120-168 MHz
        logger.info(f'{nchan_init} channels, {nchan} of which are above 168MHz')
        # nchan = (48*freqstep) * (nchan // (48*freqstep)) # multiple of 48 after average
        lib_util.check_rm('mss-shift')
        os.makedirs('mss-shift')
        logger.info(f'Phase shifting @{timeint_init}s, using {nchan} channels (from {nchan_init})')

        for MS in MSs.getListObj():
            nchan_thisms = int(np.sum(np.array(MS.getFreqs()) < 168.3e6))  # only use 120-168 MHz
            if nchan_thisms == 0:
                logger.warning(f"Skipping {MS.nameMS} (above 168.3MHz)")
                continue
            commandCurrent = MS.concretiseString(
                f'DP3 {parset_dir}/DP3-shift.parset msin=$pathMS msout=mss-shift/$nameMS.MS msin.datacolumn=CORRECTED_DATA '
                f'msin.nchan={nchan_thisms} shift.phasecenter=[{model_centre[0]}deg,{model_centre[1]}deg]')
            logCurrent = MS.concretiseString('$nameMS_initshift.log')
            s.add(cmd=commandCurrent, log=logCurrent, commandType='DP3', )
        s.run(check=True)
    MSs_shift = lib_ms.AllMSs( glob.glob('mss-shift/*.MS'), s, check_flags=False )

    with w.if_todo('apply_beam'):
        logger.info('Correcting beam: DATA -> DATA...')
        MSs_shift.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS corrbeam.invert=True', log='$nameMS_beam.log',
                      commandType='DP3')

    name_msavg = MSs_shift.getListObj()[0].getNameField() + '-combined.MS'
    if not os.path.exists('mss-shiftavg'):
        timeint = MSs_shift.getListObj()[0].getTimeInt()
        t_avg_factor = int(16/timeint)
        nchan_shift = len(MSs_shift.getFreqs())
        nchan_shiftavg = len(MSs_shift.getListObj()) / 2 # 0.5 chan/SB
        f_avg_factor = int(nchan_shift/nchan_shiftavg)
        lib_util.check_rm('mss-shiftavg')
        os.makedirs('mss-shiftavg')
        # Average all MSs to one single MS (so that we can average to less than 1 chan/SB)
        logger.info(f'Averaging {timeint}s -> {t_avg_factor*timeint}s; {nchan_shift}chan -> {nchan_shiftavg}chan')
        s.add(f'DP3 {parset_dir}/DP3-avg.parset msin=[{",".join(MSs_shift.getListStr())}] msout=mss-shiftavg/{name_msavg}'
              f' avg.timestep={t_avg_factor} avg.freqstep={f_avg_factor}', log=f'create-shiftavg.log', commandType='DP3')
        s.run(check=True)
    MSs_shiftavg = lib_ms.AllMSs(f'mss-shiftavg/{name_msavg}', s, check_flags=False )

    # Add model to MODEL_DATA
    predict_fits_model(MSs_shiftavg, fits_model)
    #####################################################################################################
    # Get mask
    if not os.path.exists(peelMask):
        logger.info('Create mask...')
        # dummy clean to get image -> mask
        lib_util.run_wsclean(s, 'wsclean-mask.log', MSs_shiftavg.getStrWsclean(), niter=0, channel_range='0 1',
                             no_reorder='', interval='0 10', name='img/mask', scale=pixscale, size=imgsize, nmiter=0)
        # create peelMask
        copy2(f'img/mask-dirty.fits', f'{peelMask}')
        lib_img.blank_image_reg(peelMask, peelReg.filename, inverse=True, blankval=0.)
        lib_img.blank_image_reg(peelMask, peelReg.filename, inverse=False, blankval=1.)

    # Self-cal cycle
    # TODO check if this already works
    if pointing_distance <= 1/3600:
        maxiter = 1
        sol_factor_t, sol_factor_f = 1, 1
    elif 1/3600 < pointing_distance <= 5:
        maxiter = 4
        sol_factor_t, sol_factor_f = 1, 1
    else:
        maxiter = 1
        sol_factor_t, sol_factor_f = 1, 1
    for c in range(maxiter):
        imagename = f'img/peel-c{c:02}'
        ##################################################################################
        # solve+apply scalarphase and fulljones
        solve_and_apply(MSs_shiftavg, f'{c:02}', sol_factor_t, sol_factor_f)
        # corrupt and subtract with above solutions
        corrupt_subtract_testimage(MSs_shiftavg, f'{c:02}')
        ##################################################################################
        # Imaging
        if c <= maxiter - 2: # SKIP last iteration imaging -> we wont use the new model anyways
            wsclean_params = {
                'scale': pixscale,
                'size': imgsize,
                'weight': 'briggs -0.6',
                'join_channels': '',
                'fit_spectral_pol': 5,
                'channels_out': 12,
                'minuv_l': uvlambdamin,
                'name': imagename,
                'no_update_model_required': '',
                'mgain': 0.85,
                'multiscale': '',
                'auto_mask': 3.0,
                'auto_threshold': 0.8,
                'baseline_averaging': 6,
                'use_wgridder': '',
            }
            with w.if_todo('imaging_c%02i' % c):
                logger.info(f'Cleaning (cycle: {c}; size: {imgsize} pix scale: {pixscale})...')
                lib_util.run_wsclean(s, f'wsclean-c{c}.log', MSs_shift.getStrWsclean(), niter=1000000,
                                     multiscale_convolution_padding=1.2, multiscale_scales='0,20,40,80,160,320',
                                     fits_mask=peelMask, **wsclean_params)
                os.system(f'cat logs/wsclean-c{c}.log | grep "background noise"')
            selfcal_model = imagename
            # predict with BLANKED model
            predict_fits_model(MSs_shiftavg, selfcal_model, stepname=f'predict_{c:02}', predict_reg=predictReg)
    # Predict last selfcal model to MSs_shift
    predict_fits_model(MSs_shift, selfcal_model, stepname=f'predict_fnal', predict_reg=predictReg)
    # Subtract the model
    corrupt_subtract_testimage(MSs_shift, 'fnal', sol_suffix=f'{c:02}') # --> SUBTRACT_DATA
    MSs_subtracted = MSs_shift
else: #CASE 2: Model and MS are aligned. Solve
    predict_fits_model(MSs, fits_model)
    solve_and_apply(MSs, field)
    corrupt_subtract_testimage(MSs, field) # --> SUBTRACTED_DATA
    MSs_subtracted = MSs
###################################################################################################################
if not os.path.exists('mss-peel'):
    timeint_init = MSs_subtracted.getListObj()[0].getTimeInt()
    avgtimeint = int(round(t_int / timeint_init))
    nchan_init = len(MSs_subtracted.getFreqs())
    nchan = np.sum(np.array(MSs_subtracted.getFreqs()) < 168.3e6)  # only use 120-168 MHz
    chan_per_sb = int(nchan_init / len(MSs_subtracted.getListStr()))
    freqstep =  int(chan_per_sb / final_chan_per_sb)
    lib_util.check_rm('mss-peel')
    os.makedirs('mss-peel')
    logger.info('Averaging in time (%is -> %is), channels: %ich -> %ich)' % (
        timeint_init, timeint_init * avgtimeint, nchan_init, nchan // freqstep))
    MSs_subtracted.run(f'DP3 {parset_dir}/DP3-shiftavg.parset msin=$pathMS msout=mss-peel/$nameMS.MS msin.datacolumn=SUBTRACTED_DATA '
            f'msin.nchan={nchan} avg.timestep={avgtimeint} avg.freqstep={freqstep} shift.phasecenter=[{phasecentre[0]}deg,{phasecentre[1]}deg] ',
            log='$nameMS_initavg.log', commandType='DP3')
    os.system("rename.ul .MS .ms mss-peel/*")
    MSs_peel = lib_ms.AllMSs(glob.glob('mss-peel/*.ms'), s, check_flags=False)
    logger.info('Correcting beam: DATA -> DATA...')
    # TODO check beam is corrected for peel direction here
    MSs_peel.run(f'DP3 {parset_dir}/DP3-beam.parset msin=$pathMS', log='$nameMS_beam.log', commandType='DP3')

MSs_peel = lib_ms.AllMSs(glob.glob('mss-peel/*.ms'), s)
# If distant pointing -> derive and apply new phase-solutions against field model
if pointing_distance > 1.0:
    with w.if_todo('init_field_model'):
        os.system(f'wget -O field.skymodel \'http://tgssadr.strw.leidenuniv.nl/cgi-bin/gsmv4.cgi?coord={phasecentre[0]},'
                  f'{phasecentre[1]}&radius={MSs_peel.getListObj()[0].getFWHM("min") / 2.}&unit=deg&deconv=y\' ')  # This assumes first MS is lowest in freq
        lsm = lsmtool.load('field.skymodel', MSs_peel.getListStr()[0])
        lsm.remove('I<0.5', applyBeam=True)
        lsm.remove(f'{peelMask} == True')  # This should remove the peel source from the field model.
        lsm.group('single', root='field')
        lsm.write('field.skymodel', clobber=True)
        lib_util.check_rm('field.skydb')
        os.system('makesourcedb outtype="blob" format="<" in=field.skymodel out=field.skydb')
        for MS in MSs_peel.getListStr():
            lib_util.check_rm(MS + '/field.skydb')
            logger.debug('Copy: field.skydb -> ' + MS)
            copy2('field.skydb', MS)
        logger.info('Predict field --> MODEL_DATA')
        MSs_peel.run(
            f'DP3 {parset_dir}/DP3-predict.parset msin=$pathMS msin.datacolumn=DATA msout.datacolumn=MODEL_DATA '
            f'pre.sourcedb=field.skydb pre.usebeammodel=True pre.sourcedb=$pathMS/field.skydb',
            log='$nameMS_pre.log', commandType='DP3')

    with w.if_todo('sol_di'): # Solve DATA vs. MODEL_DATA
        logger.info('Solving direction-independent (scalarphase)...')
        s.add(
            f'DP3 {parset_dir}/DP3-soldd.parset msin=[{",".join(MSs_peel.getListStr())}] msin.datacolumn=DATA sol.h5parm=$pathMS/di.h5 '
            f'sol.mode=scalarcomplexgain sol.nchan=4 sol.smoothnessconstraint=2.0e6 sol.solint={16 // t_int:d} '
            f'sol.uvlambdamin={uvlambdamin}',
            log=f'$nameMS_sol_di.log', commandType="DP3")
        s.run(check=True)

        lib_util.run_losoto(s, f'di', [ms + '/di.h5' for ms in MSs_peel.getListStr()], \
                            [ parset_dir+'/losoto-flag.parset',
                              parset_dir + '/losoto-plot-scalaramp.parset',
                              parset_dir + '/losoto-plot-scalarph.parset'])
        move(f'cal-di.h5', 'peel/solutions/')
        move(f'plots-di', 'peel/plots/')

    # Correct DATA -> CORRECTED_DATA
    with w.if_todo('cor_di'):
        logger.info('Direction-independent correction...')
        MSs_peel.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA cor.updateweights=False '
                      f'cor.parmdb=peel/solutions/cal-di.h5 cor.correction=phase000 cor.direction=[MODEL_DATA]', \
                      log=f'$nameMS_cor_di.log', commandType='DP3')

with w.if_todo(f'test-image-wide-c{c:02}'):
    logger.info('Cleaning wide...') # IMAGING - either on DATA or if present, CORRECTED_DATA
    lib_util.run_wsclean(s, f'wsclean-wide-c{c}.log', MSs_peel.getStrWsclean(), weight='briggs -0.5',
                         name='peel-wide',  parallel_deconvolution=1024, scale='4.0arcsec', size=5000, niter=500000,
                         join_channels='', channels_out=3, nmiter=10, fit_spectral_pol=3, minuv_l=uvlambdamin, multiscale='', multiscale_max_scales=4,
                         mgain=0.85, auto_threshold=1.0, auto_mask=4.0, no_update_model_required='', do_predict=False, local_rms='')
    os.system(f'cat logs/wsclean-wide-c{c}.log | grep "background noise"')

logger.info("Done.")
