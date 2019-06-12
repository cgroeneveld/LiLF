#!/usr/bin/env python
# -*- coding: utf-8 -*-

# initial calibration of the calibrator in circular, get and corr FR, back to linear, sol flag + effects separation

import sys, os, glob, re
import numpy as np
import lsmtool

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-3c.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_3c','parset_dir')
bl2flag = parset.get('flag','stations')
target = os.getcwd().split('/')[-1]
data_dir = '/home/fdg/lofar5/3Csurvey/%s' % target
userReg = parset.get('model','userReg')

def get_cal_dir(timestamp):
    """
    Get the proper cal directory from a timestamp
    """
    for cal_dir in sorted(glob.glob('../../cals/*/*')):
        cal_timestamps = set( [ms.split('_')[1][1:] for ms in glob.glob(cal_dir+'/cals/*MS')] )
        if timestamp in cal_timestamps:
            logger.info('Calibrator found: %s (t=%s)' % (cal_dir, timestamp))
            return cal_dir+'/cals'

    logger.error('Missing calibrator.')
    sys.exit()

##########################################################
logger.info('Cleaning...')
lib_util.check_rm('cal*h5')
lib_util.check_rm('plots*')
lib_util.check_rm('img')
os.makedirs('img')

MSs = lib_ms.AllMSs( glob.glob(data_dir+'/*MS'), s, check_flags=False)
for timestamp in set([ os.path.basename(ms).split('_')[1][1:] for ms in MSs.getListStr() ]):
    mss_toconcat = glob.glob(data_dir+'/'+target+'_t'+timestamp+'_SB*.MS')
    MS_concat = target+'_t'+timestamp+'_concat.MS'
    MS_concat_bkp = target+'_t'+timestamp+'_concat.MS-bkp'
    if os.path.exists(MS_concat_bkp): 
        logger.info('Restoring bkp data...')
        os.system('rm -r %s' % MS_concat)
        os.system('cp -r %s %s' % (MS_concat_bkp, MS_concat) )

    else:
        logger.info('Making %s...' % MS_concat)
        s.add('DPPP '+parset_dir+'/DPPP-avg.parset msin=\"'+str(mss_toconcat)+'\" msout='+MS_concat,\
            log=MS_concat+'_avg.log', commandType='DPPP')
        s.run(check=True, maxThreads=1)

        MSs = lib_ms.AllMSs( [MS_concat], s )
        
        # flag bad stations, and low-elev
        logger.info('Flagging...')
        MSs.run('DPPP '+parset_dir+'/DPPP-flag.parset msin=$pathMS msout=. ant.baseline=\"'+bl2flag+'\"', \
                    log='$nameMS_flag.log', commandType='DPPP')
        
        cal_dir = get_cal_dir(timestamp)
        h5_pa = cal_dir+'/cal-pa.h5'
        h5_amp = cal_dir+'/cal-amp.h5'
        h5_iono = cal_dir+'/cal-iono.h5'
        assert os.path.exists(h5_pa)
        assert os.path.exists(h5_amp)
        assert os.path.exists(h5_iono)
        
        # Correct fist for BP(diag)+TEC+Clock and then for beam
        
        # Apply cal sol - SB.MS:DATA -> SB.MS:CORRECTED_DATA (polalign corrected)
        logger.info('Apply solutions (pa)...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor.steps=[pa] \
                cor.pa.parmdb='+h5_pa+' cor.pa.correction=polalign', log='$nameMS_cor1.log', commandType='DPPP')
        
        # Apply cal sol - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (polalign corrected, calibrator corrected+reweight, beam corrected+reweight)
        logger.info('Apply solutions (amp/ph)...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.steps=[amp,ph] \
                cor.amp.parmdb='+h5_amp+' cor.amp.correction=amplitudeSmooth cor.amp.updateweights=True\
                cor.ph.parmdb='+h5_iono+' cor.ph.correction=phaseOrig000', log='$nameMS_cor2.log', commandType='DPPP')
        # clock?
        #MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
        #        cor.parmdb='+h5_amp+' cor.correction=amplitudeSmooth cor.updateweights=True', \
        #        log='$nameMS_cor2.log', commandType='DPPP')
        #MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
        #        cor.parmdb='+h5_iono+' cor.correction=clock000', \
        #        log='$nameMS_cor3.log', commandType='DPPP')
        #MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
        #        cor.parmdb='+h5_iono+' cor.correction=phase000', \
        #        log='$nameMS_cor4.log', commandType='DPPP')
        
        # Beam correction CORRECTED_DATA -> CORRECTED_DATA (polalign corrected, beam corrected+reweight)
        logger.info('Beam correction...')
        MSs.run('DPPP '+parset_dir+'/DPPP-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_beam.log', commandType='DPPP')

        # TODO: TEST
        # Convert to circular CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Converting to circular...')
        MSs.run('mslin2circ.py -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', log='$nameMS_circ2lin.log', commandType='python', maxThreads=10)

        # Move CORRECTED_DATA -> DATA
        logger.info('Move CORRECTED_DATA -> DATA...')
        MSs.run('taql "update $pathMS set DATA = CORRECTED_DATA"', log='$nameMS_taql.log', commandType='general')

        # bkp
        logger.info('Making backup...')
        os.system('cp -r %s %s' % (MS_concat, MS_concat_bkp) ) # do not use MS.move here as it resets the MS path to the moved one

MSs = lib_ms.AllMSs( glob.glob('*concat.MS'), s, check_flags=False )
MSs.plot_HAcov('HAcov.png')

#####################################################
# Model
logger.info('Preparing model...')
sourcedb = 'tgts.skydb'
if not os.path.exists(sourcedb):
    phasecentre = MSs.getListObj()[0].getPhaseCentre()
    fwhm = MSs.getListObj()[0].getFWHM(freq='min')
    radeg = phasecentre[0]
    decdeg = phasecentre[1]
    # get model the size of the image (radius=fwhm/2)
    os.system('wget -O tgts.skymodel "https://lcs165.lofar.eu/cgi-bin/gsmv1.cgi?coord=%f,%f&radius=%f&unit=deg"' % (radeg, decdeg, fwhm/2.)) # ASTRON
    lsm = lsmtool.load('tgts.skymodel')#, beamMS=MSs.getListObj()[0])
    lsm.remove('I<1')
    lsm.write('tgts.skymodel', clobber=True)
    os.system('makesourcedb outtype="blob" format="<" in=tgts.skymodel out=tgts.skydb')

# Predict MODEL_DATA
logger.info('Predict (DPPP)...')
MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+sourcedb, log='$nameMS_pre.log', commandType='DPPP')

###############################################################
# Selfcal
rms_noise_pre = np.inf; doamp = False
for c in range(100):

    logger.info('== Start cycle: %s ==' % c)

    logger.info('Merge source model into 1 patch...')
    lsm = lsmtool.load('tgts.skymodel')
    lsm.group('single')
    lsm.write('tgts-onepatch.skymodel', clobber=True)
    os.system('makesourcedb outtype="blob" format="<" in=tgts-onepatch.skymodel out=tgts-onepatch.skydb')
    sourcedb = 'tgts-onepatch.skydb'

    logger.info('Remove bad timestamps...')
    MSs.run( 'flagonmindata.py -f 0.5 $pathMS', log='$nameMS_flagonmindata.log', commandType='python')

    ####################################################
    # 1: Solving

    # Smooth DATA -> SMOOTHED_DATA
    logger.info('BL-based smoothing...')
    MSs.run('BLsmooth.py -r -f 0.2 -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth1.log', commandType='python')

    # First get scalar solution to fix TEC
    # solve G - group*_TC.MS:SMOOTHED_DATA
    logger.info('Solving...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-solG.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA sol.h5parm=$pathMS/calG0.h5 sol.mode=scalarcomplexgain \
            sol.antennaconstraint=[[CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA]]', \
            log='$nameMS_solG0-c'+str(c)+'.log', commandType="DPPP")
    lib_util.run_losoto(s, 'G0-c'+str(c), [ms+'/calG0.h5' for ms in MSs.getListStr()], \
                    [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-amp.parset'])

    # Correct SMOOTHED_DATA -> CORRECTED_DATA
    logger.info('Correction PH...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-cor.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA cor.parmdb=cal-G0-c'+str(c)+'.h5 cor.correction=phase000', \
            log='$nameMS_corPH-c'+str(c)+'.log', commandType='DPPP')

    # Then get phaseonly diagonal to find FR
    # solve G - group*_TC.MS:CORRECTED_DATA
    logger.info('Solving...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-solFR.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/calG1.h5 sol.mode=phaseonly \
            sol.antennaconstraint=[[CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA]]', \
            log='$nameMS_solG1-c'+str(c)+'.log', commandType="DPPP")
    lib_util.run_losoto(s, 'G1-c'+str(c), [ms+'/calG1.h5' for ms in MSs.getListStr()], \
                    [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-fr.parset'])

    # Correct DATA -> CORRECTED_DATA
#    logger.info('Correction FR...')
#    MSs.run('DPPP ' + parset_dir + '/DPPP-cor.parset msin=$pathMS cor.parmdb=cal-G1-c'+str(c)+'.h5 cor.correction=rotationmeasure000', \
#            log='$nameMS_corFR-c'+str(c)+'.log', commandType='DPPP')

#    # Smooth CORRECTED_DATA -> SMOOTHED_DATA
#    logger.info('BL-based smoothing...')
#    MSs.run('BLsmooth.py -r -f 0.2 -i CORRECTED_DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth1.log', commandType='python')

    # Get again TEC with scalr
    # solve G - group*_TC.MS:SMOOTHED_DATA
    logger.info('Solving...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-solG.parset msin=$pathMS msin.datacolumn=SMOOTHED_DATA sol.h5parm=$pathMS/calG2.h5 sol.mode=scalarcomplexgain \
            sol.antennaconstraint=[[CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA]]', \
            log='$nameMS_solG2-c'+str(c)+'.log', commandType="DPPP")
    lib_util.run_losoto(s, 'G2-c'+str(c), [ms+'/calG2.h5' for ms in MSs.getListStr()], \
                    [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-amp.parset'])

    # Correct DATA -> CORRECTED_DATA
#    logger.info('Correction FR...')
#    MSs.run('DPPP ' + parset_dir + '/DPPP-cor.parset msin=$pathMS cor.parmdb=cal-G1-c'+str(c)+'.h5 cor.correction=rotationmeasure000', \
#            log='$nameMS_corFR2-c'+str(c)+'.log', commandType='DPPP')
    logger.info('Correction PH...')
    MSs.run('DPPP ' + parset_dir + '/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb=cal-G2-c'+str(c)+'.h5 cor.correction=phase000', \
            log='$nameMS_corPH2-c'+str(c)+'.log', commandType='DPPP')

    if doamp:
        # solve G - group*_TC.MS:CORRECTED_DATA
        MSs.run('DPPP ' + parset_dir + '/DPPP-solG.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/calGa.h5 sol.mode=diagonal sol.solint=900', \
            log='$nameMS_solGa-c'+str(c)+'.log', commandType="DPPP")
        lib_util.run_losoto(s, 'Ga-c'+str(c), [ms+'/calGa.h5' for ms in MSs.getListStr()], \
                    [parset_dir+'/losoto-plot-ph.parset', parset_dir+'/losoto-plot-amp.parset', parset_dir+'/losoto-amp.parset'])

        # Correct CORRECTED_DATA -> CORRECTED_DATA
        logger.info('Correction AMP...')
        MSs.run('DPPP ' + parset_dir + '/DPPP-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.parmdb=cal-Ga-c'+str(c)+'.h5 cor.correction=amplitudeSmooth', \
            log='$nameMS_corAMP-c'+str(c)+'.log', commandType='DPPP')

    #################################################
    # 2: Cleaning
    
    logger.info('Cleaning (cycle: '+str(c)+')...')
    imagename = 'img/img-%02i' % c
    lib_util.run_wsclean(s, 'wscleanA-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=5000, scale='3arcsec', \
            weight='briggs 0.', niter=1000, no_update_model_required='', minuv_l=30, mgain=0.5, baseline_averaging=5, \
            join_channels='', fit_spectral_pol=2, channels_out=3)

    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
    im.makeMask(threshisl = 5)

    logger.info('Cleaning w/ mask (cycle: '+str(c)+')...')
    imagename = 'img/imgM-%02i' % c
    #auto_mask=5, local_rms=''
    lib_util.run_wsclean(s, 'wscleanB-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=5000, scale='3arcsec', \
            weight='briggs 0.', niter=1000000, update_model_required='', minuv_l=30, mgain=0.5, \
            fits_mask=im.maskname, multiscale='', auto_threshold=3, use_weights_as_taper='', \
            join_channels='', fit_spectral_pol=2, channels_out=3)
    os.system('cat logs/wscleanB-c'+str(c)+'.log | grep "background noise"')

    # Set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA
    logger.info('Set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA...')
    MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"', log='$nameMS_taql.log', commandType='general')

    #logger.info('Cleaning residuals (cycle: '+str(c)+')...')
    #imagename = 'img/imgR-%02i' % c
    #lib_util.run_wsclean(s, 'wscleanC-c'+str(c)+'.log', MSs.getStrWsclean(), name=imagename, size=4500, scale='4arcsec', \
    #        weight='briggs 0.', niter=1000000, update_model_required='', minuv_l=30, mgain=0.7, \
    #        join_channels='', fit_spectral_pol=2, channels_out=3)
    #os.system('cat logs/wscleanC-c'+str(c)+'.log | grep "background noise"')

    # Flag on residuals (CORRECTED_DATA)
    logger.info('Flagging residuals...')
    MSs.run('DPPP '+parset_dir+'/DPPP-flagres.parset msin=$pathMS', log='$nameMS_flagres-c'+str(c)+'.log', commandType='DPPP')

    rms_noise = lib_img.Image(imagename+'-MFS-image.fits').getNoise()
    logger.info('RMS noise: %f' % rms_noise)
    if rms_noise > 0.95*rms_noise_pre:
        if doamp: break # if already doing amp and not getting better, quit
        doamp = True
    rms_noise_pre = rms_noise

logger.info("Done.")
