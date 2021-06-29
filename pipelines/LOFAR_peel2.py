#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Peel a bright source outside of the PB from the data
# Two possible ways to run this pipeline:
# 1. You provide either a ds9 circle/polygon peelReg (default 'peel.reg') and either a sourcedb or a fits_model
#    in the [model] section. The pipeline will solve against this initial model and then continue w/ self-cal.
# 2. (not implemented yet) You provide initial solutions (e.g. from a parallel pointing), in this case the data will be corrected w/
#    these initial solutions and imaged, then it will continue with self-cal using this model.
#
# The MS files need to be in "./mss/"

import sys, os, glob, re
from shutil import copy2, copytree, move
import numpy as np
import casacore.tables as pt
import lsmtool

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-peel.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-peel.walker')

parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_peel','parset_dir')
cal_dir = parset.get('LOFAR_peel','cal_dir')
data_dir = parset.get('LOFAR_peel','data_dir')
sourcedb = parset.get('model','sourcedb')
bl2flag = parset.get('flag','stations')
#############################################################################
# TODO final beam corruption/correction?

uvlambdamin = 30
t_int = 8
freqstep = 2 # must be in [1,2,3,4,6,8,12,16,24,48]
nchunks = 20 # number of channel chunks for imaging and calibration

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

MSs = lib_ms.AllMSs( glob.glob(data_dir + '/*.MS'), s )
if not MSs.isHBA:
    logger.error('Only HBA measurement sets supported for now.')
    sys.exit(1)
try:
    MSs.print_HAcov()
except:
    logger.error('Problem with HAcov, continue anyway.')
logger.info(f'Resolution {MSs.resolution}\'\', using uvlambdamin={uvlambdamin}, averaging a factor of {freqstep} in freq and to an '
            f'integration time of t_int={t_int}')

phasecentre = MSs.getListObj()[0].getPhaseCentre()
# region must be a list of ds9 circles and polygons (other shapes can be implemented in lib_util.Rgion_helper()
# center = peelReg.get_center() # center of the extract region
center = np.rad2deg([-3.00711503, 0.21626569]) # hardcoded m87 - for a fits model, does this need to be the phase center?

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

    # Correct fist for BP(diag)+TEC+Clock and then for beam
    # Apply cal sol - SB.MS:DATA -> SB.MS:CORRECTED_DATA (polalign corrected)
    logger.info('Apply solutions (pa)...')
    MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA \
            cor.parmdb=' + h5_pa + ' cor.correction=polalign', log='$nameMS_cor1.log', commandType='DP3')

    # Apply cal sol - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (polalign corrected, calibrator corrected+reweight, beam corrected+reweight)
    logger.info('Apply solutions (amp/ph)...')
    if MSs.isLBA:
        MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.steps=[amp,ph] \
                cor.amp.parmdb=' + h5_amp + ' cor.amp.correction=amplitudeSmooth cor.amp.updateweights=True\
                cor.ph.parmdb=' + h5_iono + ' cor.ph.correction=phaseOrig000', log='$nameMS_cor2.log',
                commandType='DP3')
    elif MSs.isHBA:
        MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS cor.steps=[amp,clock] \
                cor.amp.parmdb=' + h5_amp + ' cor.amp.correction=amplitudeSmooth cor.amp.updateweights=True\
                cor.clock.parmdb=' + h5_iono + ' cor.clock.correction=clockMed000', log='$nameMS_cor2.log',
                commandType='DP3')

    # Beam correction CORRECTED_DATA -> CORRECTED_DATA (polalign corrected, beam corrected+reweight)
    logger.info('Beam correction...')
    MSs.run('DP3 ' + parset_dir + '/DP3-beam.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA '
            'msout.datacolumn=CORRECTED_DATA corrbeam.updateweights=True', log='$nameMS_beam.log', commandType='DP3')
### DONE

#################################################################################################
# Initial flagging
with w.if_todo('flag'):
    logger.info('Flagging...')
    MSs.run(f'DP3 {parset_dir}/DP3-flag.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA ant.baseline=\"{bl2flag}\" '
            f'aoflagger.strategy={parset_dir}/HBAdefaultwideband.lua uvmin.uvlambdamin={uvlambdamin}',
            log='$nameMS_flag.log', commandType='DP3')
    logger.info('Remove bad timestamps...')
    MSs.run('flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')

    logger.info('Plot weights...')
    MSs.run(f'reweight.py $pathMS -v -p -d CORRECTED_DATA -a {"CS001HBA0" if MSs.isHBA else "CS001LBA"}',
            log='$nameMS_weights.log', commandType='python')
    os.system('mv *.png self/plots')
### DONE

# Add model to MODEL_DATA
with w.if_todo('init_model'):
    assert os.path.isfile(sourcedb)
    logger.info(f'Using sourcedb {sourcedb}.')
    # copy sourcedb into each MS to prevent concurrent access from multiprocessing to the sourcedb
    sourcedb_basename = sourcedb.split('/')[-1]
    for MS in MSs.getListStr():
        lib_util.check_rm(MS + '/' + sourcedb_basename)
        logger.debug('Copy: ' + sourcedb + ' -> ' + MS)
        copy2(sourcedb, MS)

    # note: do not add MODEL_DATA or the beam is transported from DATA, while we want it without beam applied
    logger.info('Predict (DP3: %s))...' % (sourcedb_basename))
    MSs.run(f'DP3 {parset_dir}/DP3-predict.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA pre.usebeammodel=False '
            f'pre.sourcedb=$pathMS/{sourcedb_basename}', log='$nameMS_pre.log', commandType='DP3')
#####################################################################################################
### DONE
# Self-cal cycle
with w.if_todo('solve_iono'):
    # Solve cal_SB.MS: CORRECTED_DATA (only solve)
    logger.info('Solving scalarphase...')
    MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA '
            f'sol.h5parm=$pathMS/iono.h5 sol.mode=scalarcomplexgain sol.nchan=1 sol.smoothnessconstraint=3e6 '
            f'sol.uvlambdamin={uvlambdamin}' , log=f'$nameMS_sol_iono.log', commandType="DP3")

    lib_util.run_losoto(s, f'iono', [ms+'/iono.h5' for ms in MSs.getListStr()], \
                        [#parset_dir+'/losoto-flag.parset',
                         parset_dir+'/losoto-plot-scalaramp.parset',
                         parset_dir+'/losoto-plot-scalarph.parset'])

    move(f'cal-iono.h5', 'peel/solutions/')
    move(f'plots-iono', 'peel/plots/')

with w.if_todo('cor_iono'):
    logger.info('Scalarphase correction...')
    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA msout.datacolumn=CORRECTED_DATA2 cor.updateweights=False '
            f'cor.parmdb=peel/solutions/cal-iono.h5 cor.correction=phase000', \
            log=f'$nameMS_cor_iono.log', commandType='DP3')
    ### DONE
with w.if_todo('solve_gain'):
    logger.info('Solving gain...')
    MSs.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA2 sol.h5parm=$pathMS/gain.h5 '
            f'sol.mode=rotation+diagonal sol.smoothnessconstraint=4e6 sol.nchan=4 sol.solint={300 // t_int} '
            f'sol.uvlambdamin={uvlambdamin}', log=f'$nameMS_sol_gain.log', commandType="DP3")

    lib_util.run_losoto(s, f'gain', [ms + '/gain.h5' for ms in MSs.getListStr()], [parset_dir + '/losoto-diag.parset',
                                                                                   parset_dir + '/losoto-plot-rot.parset']) # maybe remove norm for subtract?

    move(f'cal-gain.h5', 'peel/solutions/')
    move(f'plots-gain', 'peel/plots/')

# Correct gain amp and ph CORRECTED_DATA -> CORRECTED_DATA
# with w.if_todo('cor_gain'):
#     logger.info('Gain correction...')
#     MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.steps=[amp,phase] \
#             cor.amp.parmdb=peel/solutions/cal-gain.h5 cor.amp.correction=amplitude000 '
#             f'cor.parmdb=peel/solutions/cal-gain.h5 cor.phase.correction=phase000',
#             log=f'$nameMS_cor_gain.log',
#             commandType='DP3')

# with w.if_todo('solve_fulljones%02i' % c):
#     logger.info('Solving full-Jones...')
#     MSs_peel.run(f'DP3 {parset_dir}/DP3-soldd.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA '
#                  f'sol.h5parm=$pathMS/fulljones.h5 sol.mode=fulljones sol.nchan=1 sol.solint={300//t_int} '
#                  f'sol.smoothnessconstraint=4e6 sol.uvlambdamin={uvlambdamin}',
#                  log=f'$nameMS_sol_fulljones-c{c}.log', commandType="DP3")
#
#     lib_util.run_losoto(s, f'fulljones-c{c:02}', [ms + '/fulljones.h5' for ms in MSs_peel.getListStr()], \
#                         [#parset_dir + '/losoto-norm.parset',
#                          parset_dir + '/losoto-plot-fulljones.parset'])
#
#     move(f'cal-fulljones-c{c:02}.h5', 'peel/solutions/')
#     move(f'plots-fulljones-c{c:02}', 'peel/plots/')
#
# # Correct gain amp and ph CORRECTED_DATA -> CORRECTED_DATA
# with w.if_todo('cor_fulljones_c%02i' % c):
#     logger.info('Full-Jones correction...')
#     MSs_peel.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS cor.correction=fulljones '
#             f'cor.parmdb=peel/solutions/cal-fulljones-c{c:02}.h5 cor.soltab=\[amplitude000,phase000\]',
#             log=f'$nameMS_cor_gain-c{c:02}.log', commandType='DP3')

###################################################################################################################

with w.if_todo(f'subtract'):
    logger.info('Scalarphase corruption...')
    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA '
                 f'cor.updateweights=False cor.parmdb=peel/solutions/cal-iono.h5 cor.correction=phase000 cor.invert=False', \
                 log=f'$nameMS_corrup_iono.log', commandType='DP3')
    logger.info('Gain corruption...')
    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=MODEL_DATA '
                 f'msout.datacolumn=MODEL_DATA cor.steps=[amp,phase,rot] cor.invert=False '
                 f'cor.amp.parmdb=peel/solutions/cal-gain.h5 cor.amp.correction=amplitude000 '
                 f'cor.parmdb=peel/solutions/cal-gain.h5 cor.phase.correction=phase000 '
                 f'cor.rot.parmdb=peel/solutions/cal-gain.h5 cor.rot.correction=rotation000',
                 log=f'$nameMS_cor_gain.log',
                 commandType='DP3')
    # logger.info('Full-Jones corruption...')
    # MSs_peel.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRUPTED_MODEL_DATA '
    #              f'msout.datacolumn=CORRUPTED_MODEL_DATA cor.correction=fulljones '
    #              f'cor.parmdb=peel/solutions/cal-fulljones-c{c:02}.h5 cor.soltab=\[phase000,amplitude000\] '
    #              f'cor.invert=False', log=f'$nameMS_corrup_gain-c{c:02}.log', commandType='DP3')
    # Set MODEL_DATA = 0 where data are flagged, then unflag everything

    MSs.addcol('SUBTRACTED_DATA', 'CORRECTED_DATA')
    logger.info('SET SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA') # CORRECTED_DATA is corrected for calibrator + beam (phase center)
    MSs.run('taql "UPDATE $pathMS SET SUBTRACTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql_subtract.log',
                 commandType='general')

with w.if_todo('die_correct'):
    logger.info('Scalarphase correction...')
    logger.info('Correcting scalarphase: SUBTRACTED_DATA -> SUBTRACTED_DATA...')
    MSs.run(f'DP3 {parset_dir}/DP3-cor.parset msin=$pathMS msin.datacolumn=SUBTRACTED_DATA msout.datacolumn=SUBTRACTED_DATA cor.updateweights=False '
                 f'cor.parmdb=peel/solutions/cal-iono.h5 cor.correction=phase000', \
                 log=f'$nameMS_cor_iono.log', commandType='DP3')


if not os.path.exists('mss-peel'):
    timeint_init = MSs.getListObj()[0].getTimeInt()
    avgtimeint = int(round(t_int/timeint_init))  # to 16 seconds
    nchan_init = len(MSs.getFreqs())
    nchan = np.sum(np.array(MSs.getFreqs()) < 168.3e6) # only use 120-168 MHz
    logger.info(f'{nchan_init} channels, {nchan} of which are above 168MHz')
    nchan = (48*freqstep) * (nchan // (48*freqstep)) # multiple of 48 after average
    lib_util.check_rm('mss-peel')
    os.makedirs('mss-peel')
    logger.info('Averaging in time (%is -> %is), channels: %ich -> %ich)' % (timeint_init,timeint_init*avgtimeint,nchan_init,nchan//freqstep))
    MSs.run(f'DP3 {parset_dir}/DP3-avg.parset msin=$pathMS msout=mss-peel/$nameMS.MS msin.datacolumn=SUBTRACTED_DATA '
            f'msin.nchan={nchan} avg.timestep={avgtimeint} avg.freqstep={freqstep}',
            log='$nameMS_avg.log', commandType='DP3')

with w.if_todo(f'clean-wide'):
    logger.info('Test empty wide... (SUBTRACTED_DATA)')
    MSs_peel = lib_ms.AllMSs(glob.glob('mss-peel/*.MS'), s)
    lib_util.run_wsclean(s, f'wsclean-peel.log', MSs_peel.getStrWsclean(), weight='briggs -0.5', data_column='DATA',
                         name='peel-wide-empty', scale='4.0arcsec', size=7000, niter=0, nmiter=0, minuv_l=uvlambdamin)
logger.info("Done.")
