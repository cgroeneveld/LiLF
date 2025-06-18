#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, glob
import numpy as np
import lsmtool
from astropy.table import Table as astrotab

########################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log
logger_obj = lib_log.Logger('pipeline-3c.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-3c.walker')

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_3c', 'parset_dir')
skydb_demix = parset.get('LOFAR_demix','demix_model')
bl2flag = parset.get('flag', 'stations')

target = os.getcwd().split('/')[-1]
data_dir = '/home/fdg/lofar5/3Csurvey/%s' % target
extended_targets = ['3c223','3c231','3c236','3c264','3c274','3c284','3c285','3c293','3c296','3c31','3c310',
                    '3c326','3c33','3c35','3c382','3c386','3c442a','3c449','3c454.3','3c465','3c84']
very_extended_targets = ['3c138','da240']

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
with w.if_todo('setup'):
    logger.info('Cleaning...')
    lib_util.check_rm('cal*h5')
    lib_util.check_rm('plots*')
    lib_util.check_rm('peel*')
    lib_util.check_rm('img')
    os.makedirs('img')

    MSs = lib_ms.AllMSs( glob.glob(data_dir+'/*MS'), s, check_flags=False)
    for timestamp in set([ os.path.basename(ms).split('_')[1][1:] for ms in MSs.getListStr() ]):
        mss_toconcat = sorted(glob.glob(data_dir+'/'+target+'_t'+timestamp+'_SB*.MS'))
        MS_concat = target+'_t'+timestamp+'_concat.MS'
        MS_concat_bkp = target+'_t'+timestamp+'_concat.MS-bkp'
    
        if os.path.exists(MS_concat_bkp): 
            logger.info('Restoring bkp data: %s...' % MS_concat_bkp)
            lib_util.check_rm(MS_concat)
            os.system('cp -r %s %s' % (MS_concat_bkp, MS_concat) )
    
        else:
            logger.info('Making %s...' % MS_concat)
            s.add('DP3 '+parset_dir+'/DP3-avg.parset msin=\"'+str(mss_toconcat)+'\" msout='+MS_concat,\
                log=MS_concat+'_avg.log', commandType='DP3')
            s.run(check=True, maxProcs=1)
    
            MSs = lib_ms.AllMSs( [MS_concat], s )
            
            # flag bad stations, and low-elev
            logger.info('Flagging...')
            MSs.run('DP3 '+parset_dir+'/DP3-flag.parset msin=$pathMS msout=. \
                     aoflagger.strategy='+parset_dir+'/LBAdefaultwideband.rfis ant.baseline=\"'+bl2flag+'\"',
                     log='$nameMS_flag.log', commandType='DP3')
            
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
            MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS \
                    cor.parmdb='+h5_pa+' cor.correction=polalign', log='$nameMS_cor1.log', commandType='DP3')
            
            # Apply cal sol - SB.MS:CORRECTED_DATA -> SB.MS:CORRECTED_DATA (polalign corrected, calibrator corrected+reweight, beam corrected+reweight)
            logger.info('Apply solutions (amp/ph)...')
            MSs.run('DP3 '+parset_dir+'/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA cor.steps=[amp,ph] \
                    cor.amp.parmdb='+h5_amp+' cor.amp.correction=amplitudeSmooth cor.amp.updateweights=True\
                    cor.ph.parmdb='+h5_iono+' cor.ph.correction=phaseOrig000', log='$nameMS_cor2.log', commandType='DP3')
            
            # Beam correction CORRECTED_DATA -> CORRECTED_DATA (polalign corrected, beam corrected+reweight)
            logger.info('Beam correction...')
            MSs.run('DP3 '+parset_dir+'/DP3-beam.parset msin=$pathMS corrbeam.updateweights=True', log='$nameMS_beam.log', commandType='DP3')
    
            # Convert to circular CORRECTED_DATA -> CORRECTED_DATA
            #logger.info('Converting to circular...')
            #MSs.run('mslin2circ.py -i $pathMS:CORRECTED_DATA -o $pathMS:CORRECTED_DATA', log='$nameMS_circ2lin.log', commandType='python', maxProcs=10)
    
            # Move CORRECTED_DATA -> DATA
            logger.info('Move CORRECTED_DATA -> DATA...')
            MSs.run('taql "update $pathMS set DATA = CORRECTED_DATA"',
                    log='$nameMS_taql.log', commandType='general')

            # bkp
            logger.info('Making backup...')
            os.system('cp -r %s %s' % (MS_concat, MS_concat_bkp) ) # do not use MS.move here as it resets the MS path to the moved one

    MSs_orig = lib_ms.AllMSs( glob.glob('*concat.MS'), s, check_flags=False )

    # Demix
    for ateam in ['VirA', 'TauA', 'CygA', 'CasA']:  #, '3C338', '3C380']:
        sep = MSs_orig.getListObj()[0].distBrightSource(ateam)
        logger.info('%s - sep: %.0f deg' % (ateam, sep))
        if sep > 2 and sep < 25 and (ateam != 'CasA' and ateam != 'CygA'):
            logger.warning('Demix of %s (sep: %.1f deg)' % (ateam, sep))
            for MS in MSs_orig.getListStr():
                lib_util.check_rm(MS + '/' + os.path.basename(skydb_demix))
                os.system('cp -r ' + skydb_demix + ' ' + MS + '/' + os.path.basename(skydb_demix))

            # TODO make a single patch for source skymodel and use that in the demix?
            logger.info('Demixing...')
            MSs_orig.run('DP3 ' + parset_dir + '/DP3-demix.parset msin=$pathMS msout=$pathMS demixer.skymodel=$pathMS/' + os.path.basename(skydb_demix) +
                ' demixer.instrumentmodel=$pathMS/instrument_demix demixer.subtractsources=[' + ateam + ']',
                log='$nameMS_demix.log', commandType='DP3')

    # Phase up stations DATA -> DATA
    lib_util.check_rm('*MS-phaseup')
    logger.info('Phase up superterp DATA -> DATA...')
    MSs_orig.run('DP3 '+parset_dir+'/DP3-phaseup.parset msin=$pathMS msout=$pathMS-phaseup',
                 log='$nameMS_phaseup.log', commandType='DP3')
    os.system('rm -r *concat.MS')
### DONE

MSs = lib_ms.AllMSs( glob.glob('*concat.MS-phaseup'), s, check_flags=False, check_sun=True )
MSs.print_HAcov('plotHAelev.png')
MSs.getListObj()[0].makeBeamReg('beam02.reg', freq='mid', pb_cut=0.2)
beam02Reg = 'beam02.reg'
MSs.getListObj()[0].makeBeamReg('beam07.reg', freq='mid', pb_cut=0.7)
beam07reg = 'beam07.reg'
region = '%s/regions/%s.reg' % (parset_dir, target)
if not os.path.exists(region): region = None

#####################################################
# Model
with w.if_todo('predict'):
    logger.info('Preparing model...')
    sourcedb = 'tgts.skydb'
    if not os.path.exists(sourcedb):
        phasecentre = MSs.getListObj()[0].getPhaseCentre()
        fwhm = MSs.getListObj()[0].getFWHM(freq='min')
        radeg = phasecentre[0]
        decdeg = phasecentre[1]
        # get model the size of the image (radius=fwhm/2)
        os.system('wget -O tgts.skymodel "https://lcs165.lofar.eu/cgi-bin/gsmv1.cgi?coord=%f,%f&radius=%f&unit=deg"' % (radeg, decdeg, fwhm)) # ASTRON
        lsm = lsmtool.load('tgts.skymodel')#, beamMS=MSs.getListObj()[0])
        lsm.remove('I<0.5')
        lsm.write('tgts.skymodel', clobber=True)
        os.system('makesourcedb outtype="blob" format="<" in=tgts.skymodel out=tgts.skydb')
    
    # Predict MODEL_DATA
    logger.info('Predict (DP3)...')
    MSs.run('DP3 '+parset_dir+'/DP3-predict.parset msin=$pathMS pre.sourcedb='+sourcedb, log='$nameMS_pre.log', commandType='DP3')
    
    # Smooth DATA -> DATA
    logger.info('BL-based smoothing...')
    MSs.run('BLsmooth.py -r -s 0.8 -i DATA -o DATA $pathMS', log='$nameMS_smooth1.log', commandType='python')
### DONE

###############################################################
# Selfcal
rms_noise_pre = np.inf; mm_ratio_pre = 0; doamp = False
solint_ph = lib_util.Sol_iterator([10,3,1])
solint_amp = lib_util.Sol_iterator([200,100,50])
for c in range(100):
    logger.info('== Start cycle: %s ==' % c)

    ####################################################
    # 1: Solving

    with w.if_todo('calib-ph-c%02i' % c):
        # solve G - group*_TC.MS:CORRECTED_DATA
        logger.info('Solving fast...')
        solint = next(solint_ph)
        MSs.run('DP3 ' + parset_dir + '/DP3-solG.parset msin=$pathMS msin.datacolumn=DATA sol.h5parm=$pathMS/calGp.h5 sol.mode=scalarcomplexgain \
                sol.solint='+str(solint)+' sol.smoothnessconstraint=5e6',
                log='$nameMS_solGp-c'+str(c)+'.log', commandType="DP3")
        lib_util.run_losoto(s, 'Gp-c'+str(c), [ms+'/calGp.h5' for ms in MSs.getListStr()],
                        [parset_dir+'/losoto-clip-large.parset', parset_dir+'/losoto-plot2d.parset', parset_dir+'/losoto-plot.parset'])
    
        # Correct DATA -> CORRECTED_DATA
        logger.info('Correction PH...')
        MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=DATA cor.parmdb=cal-Gp-c'+str(c)+'.h5 cor.correction=phase000', \
                log='$nameMS_corPH-c'+str(c)+'.log', commandType='DP3')
    ### DONE

    if doamp:
        with w.if_todo('calib-amp-c%02i' % c):
            # solve G - group*_TC.MS:CORRECTED_DATA

            logger.info('Solving slow...')
            #sol.antennaconstraint=[[CSsuperLBA,CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA,RS106LBA,RS205LBA,RS208LBA,RS210LBA,RS305LBA,RS306LBA,RS307LBA,RS310LBA,RS406LBA,RS407LBA,RS409LBA,RS503LBA,RS508LBA,RS509LBA]] \
            solint = next(solint_amp)
            MSs.run('DP3 ' + parset_dir + '/DP3-solG.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA sol.h5parm=$pathMS/calGa.h5 sol.mode=fulljones \
                    sol.solint='+str(solint)+' sol.smoothnessconstraint=2e6',
                    log='$nameMS_solGa-c'+str(c)+'.log', commandType="DP3")
            lib_util.run_losoto(s, 'Ga-c'+str(c), [ms+'/calGa.h5' for ms in MSs.getListStr()],
                        [parset_dir+'/losoto-clip.parset', parset_dir+'/losoto-plot2d.parset', parset_dir+'/losoto-plot2d-pol.parset', parset_dir+'/losoto-plot-pol.parset'])
                        #, parset_dir+'/losoto-ampnorm.parset'])
    
            # Correct CORRECTED_DATA -> CORRECTED_DATA
            logger.info('Correction slow AMP+PH...')
            MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS msin.datacolumn=CORRECTED_DATA \
                    cor.parmdb=cal-Ga-c'+str(c)+'.h5 cor.correction=fulljones cor.soltab=\[amplitude000,phase000\]',
                    log='$nameMS_corAMPPHslow-c'+str(c)+'.log', commandType='DP3')
        ### DONE

    #################################################
    # 2: Sub field + peel
    if c == 100:
        with w.if_todo('sub-field'):

            # Low res image
            logger.info('Cleaning wide 1...')
            imagename = 'img/img-wide'
            lib_util.run_wsclean(s, 'wsclean-wide.log', MSs.getStrWsclean(), name=imagename,
                             parallel_gridding=4, baseline_averaging='', size=2500,
                             scale='10arcsec', weight='briggs -0.7',
                             taper_gaussian='30arcsec',
                             niter=1000000, no_update_model_required='', minuv_l=30, mgain=0.75, nmiter=0,
                             auto_threshold=5,  local_rms='', local_rms_method='rms-with-min',
                             join_channels='', fit_spectral_pol=2, channels_out=2)
            os.system('cat logs/wsclean-wide.log | grep "background noise"')

            # makemask
            im = lib_img.Image(imagename + '-MFS-image.fits', userReg=region)
            im.makeMask(threshpix=5, rmsbox=(50, 5))
            maskfits = imagename + '-mask.fits'

            logger.info('Cleaning wide 2...')
            imagenameM = 'img/img-wideM'
            lib_util.run_wsclean(s, 'wsclean-wide.log', MSs.getStrWsclean(), name=imagenameM, do_predict=True,
                                 parallel_gridding=4, baseline_averaging='', size=2500, reuse_psf=imagename, reuse_dirty=imagename,
                                 scale='10arcsec', weight='briggs -0.7',
                                 taper_gaussian='30arcsec', fits_mask=maskfits,
                                 niter=1000000, no_update_model_required='', minuv_l=30, mgain=0.75, nmiter=0,
                                 auto_threshold=1, auto_mask=3, local_rms='', local_rms_method='rms-with-min',
                                 join_channels='', fit_spectral_pol=2, channels_out=2)

            # subtract everything
            logger.info('Subtract model: CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA...')
            MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"',
                    log='$nameMS_taql.log', commandType='general')
        # DONE

        # load skymodel
        imagename = 'img/img-wideM'
        full_image = lib_img.Image(imagename+'-MFS-image.fits', userReg=region)
        mask_ddcal = full_image.imagename.replace('.fits', '_mask-ddcal.fits')  # this is used to find calibrators
        full_image.makeMask(threshpix=5, atrous_do=False, maskname=mask_ddcal, write_srl=True, write_ds9=True)
        cal = astrotab.read(mask_ddcal.replace('fits', 'cat.fits'), format='fits')
        cal = cal[np.where(cal['Total_flux'] > 3)]
        cal.sort('Total_flux')
        cal.reverse()

        # cycle on sources to peel
        phasecentre = MSs.getListObj()[0].getPhaseCentre()
        for peelsou in cal:
            name = str(peelsou['Source_id'])
            # Skip if source is close to phase centre
            dist = lib_util.distanceOnSphere(phasecentre[0], phasecentre[1], peelsou['RA'], peelsou['DEC'])
            if dist < 0.5:
                continue

            logger.info('Peeling %s (%.1f Jy)' % (name, peelsou['Total_flux']))

            with w.if_todo('peel-%s' % name):
                lib_util.check_rm('peel-%s' % name)
                os.system('mkdir peel-%s' % name)

                # make a region
                from pyregion.parser_helper import Shape
                import pyregion
                peel_region_file = 'peel-'+name+'/'+name+'.reg'
                sh = Shape('circle', None)
                sh.coord_format = 'fk5'
                sh.coord_list = [peelsou['RA'], peelsou['DEC'], 0.1]  # ra, dec, diam
                sh.coord_format = 'fk5'
                sh.attr = ([], {'width': '2', 'point': 'cross', 'font': '"helvetica 16 normal roman"'})
                sh.comment = 'color=red text="%s"' % (name+'.reg')
                regions = pyregion.ShapeList([sh])
                lib_util.check_rm(peel_region_file)
                regions.write(peel_region_file)

                # copy and blank models
                logger.info('Peel - Cleanup model images...')
                os.system('cp '+imagenameM+'*model.fits peel-'+name)
                imagename_peel = 'peel-'+name+'/'+imagenameM.split('/')[-1]
                for model_file in glob.glob(imagename_peel + '*model.fits'):
                    lib_img.blank_image_reg(model_file, peel_region_file, blankval=0., inverse=True)

                # predict the source to peel
                logger.info('Peel - Predict init...')
                s.add('wsclean -predict -name ' + imagename_peel + ' -j ' + str(s.max_cpucores) + ' -channels-out 2 \
                      -reorder -parallel-reordering 4 ' + MSs.getStrWsclean(),
                      log='wsclean-pre.log', commandType='wsclean')
                s.run(check=True)

                # add the source to peel back
                logger.info('Peel - add model: CORRECTED_DATA = CORRECTED_DATA + MODEL_DATA...')
                MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA + MODEL_DATA"',
                        log='$nameMS_taql.log', commandType='general')

                # phaseshift + avg
                logger.info('Peel - Phaseshift+avg...')
                lib_util.check_rm('mss-dir')
                os.makedirs('mss-dir')
                MSs.run('DP3 '+parset_dir+'/DP3-shiftavg.parset msin=$pathMS msout=mss-dir/$nameMS.MS \
                    msin.datacolumn=CORRECTED_DATA msout.datacolumn=DATA \
                    avg.timestep=8 avg.freqstep=16 \
                    shift.phasecenter=\['+str(peelsou['RA'])+'deg,'+str(peelsou['DEC'])+'deg\]', \
                    log='$nameMS_shift.log', commandType='DP3')
                MSs_shift = lib_ms.AllMSs(glob.glob('mss-dir/*.MS'), s, check_flags=False, check_sun=True)

                # image
                logger.info('Peel - Image...')
                imagename_peel = 'peel-%s/img_%s' % (name, name)
                lib_util.run_wsclean(s, 'wsclean-c%02i-peel.log' % c, MSs_shift.getStrWsclean(),
                                     do_predict=True, name=imagename_peel, size=512,
                                     parallel_gridding=4, baseline_averaging='', scale='2.5arcsec',
                                     niter=100000, no_update_model_required='', minuv_l=30, mgain=0.4, nmiter=0,
                                     auto_threshold=5, local_rms='', local_rms_method='rms-with-min',
                                     join_channels='', fit_spectral_pol=2, channels_out=2)

                # calibrate
                logger.info('Peel - Calibrate...')
                MSs_shift.run('DP3 ' + parset_dir + '/DP3-solG.parset msin=$pathMS msin.datacolumn=DATA \
                        sol.h5parm=$pathMS/calGp.h5 sol.mode=scalarcomplexgain \
                        sol.solint=10 sol.smoothnessconstraint=5e6',
                        log='$nameMS_solGp-peel.log', commandType="DP3")
                lib_util.run_losoto(s, 'Gp-peel_%s' % name,
                                    [ms + '/calGp.h5' for ms in MSs_shift.getListStr()],
                                    [parset_dir + '/losoto-plot2d.parset', parset_dir + '/losoto-plot.parset'],
                                    plots_dir='peel-%s' % name)

                # predict in MSs
                logger.info('Peel - Predict final...')
                for model_file in glob.glob(imagename_peel + '*model.fits'):
                    lib_img.blank_image_reg(model_file, peel_region_file, blankval=0., inverse=True)
                s.add('wsclean -predict -name ' + imagename_peel + ' -j ' + str(s.max_cpucores) + ' -channels-out 2 \
                      -reorder -parallel-reordering 4 ' + MSs.getStrWsclean(),
                      log='wsclean-pre.log', commandType='wsclean')
                s.run(check=True)

                # corrupt
                MSs.run('DP3 ' + parset_dir + '/DP3-cor.parset msin=$pathMS \
                        msin.datacolumn=MODEL_DATA msout.datacolumn=MODEL_DATA \
                        cor.invert=False cor.parmdb=cal-Gp-peel_%s.h5 cor.correction=phase000' % name,
                        log='$nameMS_corrupt.log', commandType='DP3')

                # subtract
                logger.info('Subtract model: CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA...')
                MSs.run('taql "update $pathMS set CORRECTED_DATA = CORRECTED_DATA - MODEL_DATA"',
                    log='$nameMS_taql.log', commandType='general')
            # DONE

        with w.if_todo('reprepare dataset'):
            # blank models
            logger.info('Cleanup model images...')
            for model_file in glob.glob(imagename+'*model.fits'):
                lib_img.blank_image_reg(model_file, beam07reg, blankval=0., inverse=True)

            # ft models
            s.add('wsclean -predict -name '+imagename+' -j '+str(s.max_cpucores)+' -channels-out 2 \
                  -reorder -parallel-reordering 4 '+MSs.getStrWsclean(),
                  log='wsclean-pre.log', commandType='wsclean')
            s.run(check=True)

            # prepare new data
            logger.info('Subtract model: DATA = CORRECTED_DATA + MODEL_DATA...')
            MSs.run('taql "update $pathMS set DATA = CORRECTED_DATA + MODEL_DATA"',
                    log='$nameMS_taql.log', commandType='general')
        # DONE

    #################################################
    # 3: Cleaning
    imagename = 'img/img-%02i' % c

    # TODO: average down in time/freq?

    with w.if_todo('image-c%02i' % c):
        # special for extended sources:
        if target in very_extended_targets:
            kwargs1 = {'weight': 'briggs -0.5', 'taper_gaussian': '75arcsec', 'multiscale': '', 'multiscale_scale_bias':0.5, 'multiscale_scales':'0,30,60,120,340'}
            kwargs2 = {'weight': 'briggs -0.5', 'taper_gaussian': '75arcsec', 'multiscale_scales': '0,30,60,120,340'}
        elif target in extended_targets:
            kwargs1 = {'weight': 'briggs -0.7', 'taper_gaussian': '25arcsec'}
            kwargs2 = {'weight': 'briggs -0.7', 'taper_gaussian': '25arcsec', 'multiscale_scales': '0,15,30,60,120,240'}
        else:
            kwargs1 = {'weight': 'briggs -0.8'}
            kwargs2 = {'weight': 'briggs -0.8', 'multiscale_scales': '0,10,20,40,80,160'}

        if c < 1:
            kwargs1['size'] = 4000
            kwargs2['size'] = 4000
        else:
            kwargs1['size'] = 2000
            kwargs2['size'] = 2000

        # if next is a "cont" then I need the do_predict
        logger.info('Cleaning shallow (cycle: '+str(c)+')...')
        lib_util.run_wsclean(s, 'wsclean-c%02i.log' % c, MSs.getStrWsclean(), do_predict=True, name=imagename,
                parallel_gridding=4, baseline_averaging='', scale='2.5arcsec',
                niter=1000, no_update_model_required='', minuv_l=30, mgain=0.4, nmiter=0,
                auto_threshold=5, local_rms='', local_rms_method='rms-with-min',
                join_channels='', fit_spectral_pol=2, channels_out=2, **kwargs1)

        # check if hand-made mask is available
        im = lib_img.Image(imagename+'-MFS-image.fits')
        im.makeMask( threshpix=5, rmsbox=(50,5), atrous_do=True )
        maskfits = imagename+'-mask.fits'
        if region is not None:
            lib_img.blank_image_reg(maskfits, beam02Reg, blankval = 0.)
            lib_img.blank_image_reg(maskfits, region, blankval = 1.)

        logger.info('Cleaning full (cycle: '+str(c)+')...')
        lib_util.run_wsclean(s, 'wsclean-c%02i.log' % c, MSs.getStrWsclean(), do_predict=True, cont=True, name=imagename,
                parallel_gridding=4, scale='2.5arcsec',
                niter=1000000, no_update_model_required='', minuv_l=30, mgain=0.4, nmiter=0,
                auto_threshold=0.5, auto_mask=2., local_rms='', local_rms_method='rms-with-min', fits_mask=maskfits,
                multiscale='', multiscale_scale_bias=0.8,
                join_channels='', fit_spectral_pol=2, channels_out=2, **kwargs2)
        os.system('cat logs/wsclean-c%02i.log | grep "background noise"' % c)
    ### DONE

    im = lib_img.Image(imagename+'-MFS-image.fits')
    im.makeMask( threshpix=5, rmsbox=(500,30), atrous_do=False )
    rms_noise = im.getNoise(); mm_ratio = im.getMaxMinRatio()
    logger.info('RMS noise: %f - MM ratio: %f' % (rms_noise, mm_ratio))
    if doamp and rms_noise > 0.99*rms_noise_pre and mm_ratio < 1.01*mm_ratio_pre and c > 6:
        break  # if already doing amp and not getting better, quit
    if rms_noise > 0.95*rms_noise_pre and mm_ratio < 1.05*mm_ratio_pre:
        doamp = True
    rms_noise_pre = rms_noise; mm_ratio_pre = mm_ratio

# Low res image
logger.info('Cleaning low-res...')
imagename = 'img/img-low'
lib_util.run_wsclean(s, 'wsclean-lr.log', MSs.getStrWsclean(), name=imagename, save_source_list='',
        parallel_gridding=4, size=500, scale='10arcsec', weight='briggs -0.7', taper_gaussian='60arcsec',
        niter=1000000, no_update_model_required='', minuv_l=30, mgain=0.75, nmiter=0,
        auto_threshold=0.5, auto_mask=1, local_rms='',
        multiscale='', multiscale_scale_bias=0.8, multiscale_scales='0,10,20,40,80,160',
        join_channels='', fit_spectral_pol=2, channels_out=2)
os.system('cat logs/wsclean-lr.log | grep "background noise"')

logger.info("Done.")
