#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline for direction dependent calibration

import sys, os, glob, re
import numpy as np
import pyrap.tables as pt
import lsmtool

#######################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd
lib_log.set_logger('pipeline-dd.logger')
logger = lib_log.logger
s = lib_util.Scheduler(dry = False)

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('dd','parset_dir')
maxniter = parset.getint('dd','maxniter')
userReg = parset.get('model','userReg')

####################################################
MSs = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s )

# make beam
phasecentre = MSs.getListObj()[0].getPhaseCentre()
MSs.getListObj()[0].makeBeamReg('self/beam.reg', to_null=True) # SPARSE: go to 12 deg, first null - OUTER: go to 7 deg, first null
beamReg = 'self/beam.reg'

##########################
logger.info('Cleaning...')
lib_util.check_rm('ddcal')
os.makedirs('ddcal/masks')
os.makedirs('ddcal/plots')
os.makedirs('ddcal/images')
os.makedirs('ddcal/skymodels')

def clean(c, MSs, size=2.):
    """
    c = cycle/name
    mss = list of mss to clean
    size = in deg of the image
    """
    # set pixscale and imsize
    pixscale = MSs[0].getResolution()/3.
    imsize = int(size/(pixscale/3600.))

    if imsize < 512:
        imsize = 512

    if imsize % 2 == 1: imsize += 1 # make even

    logger.debug('Image size: '+str(imsize)+' - Pixel scale: '+str(pixscale))

    # clean 1
    logger.info('Cleaning ('+str(c)+')...')
    imagename = 'img/ddcal-'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size '+str(imsize)+' '+str(imsize)+' \
            -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale '+str(pixscale)+'arcsec -weight briggs 0. -niter 100000 -no-update-model-required -mgain 0.9 -pol I \
            -joinchannels -fit-spectral-pol 2 -channelsout 10 \
            -auto-threshold 20 -minuv-l 30 '+' '.MSs.getStrWsclean(), \
            log='wsclean-c'+str(c)+'.log', commandType='wsclean', processors='max')
    s.run(check=True)

    # make mask
    im = lib_img.Image(imagename+'-MFS-image.fits', userReg=userReg)
    im.makeMask(threshisl = 3)

    # clean 2
    #-multiscale -multiscale-scale-bias 0.5 \
    #-auto-mask 3 -rms-background-window 40 -rms-background-method rms-with-min \
    logger.info('Cleaning w/ mask ('+str(c)+')...')
    imagename = 'img/ddcalM-'+str(c)
    s.add('wsclean -reorder -name ' + imagename + ' -size '+str(imsize)+' '+str(imsize)+' \
            -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
            -scale '+str(pixscale)+'arcsec -weight briggs 0. -niter 1000000 -no-update-model-required -mgain 0.8 -pol I \
            -joinchannels -fit-spectral-pol 2 -channelsout 10 \
            -auto-threshold 0.1 -save-source-list -minuv-l 30 -fits-mask '+im.maskname+' '+' '.MSs.getStrWsclean(), \
            log='wscleanM-c'+str(c)+'.log', commandType='wsclean', processors='max')
    s.run(check=True)
    os.system('cat logs/wscleanM-c'+str(c)+'.log | grep "background noise"')

    return imagename


############################################################
logger.info('Copy data...')
if not os.path.exists('mss-dd'):
    os.makedirs('mss-dd')
    MSs.run('DPPP '+parset_dir+'/DPPP-avg.parset msin=$pathMS msout=mss-dd/$nameMS.MS msin.datacolumn=CORRECTED_DATA avg.freqstep=1 avg.timestep=1', \
                log='$nameMS_avg.log', commandType='DPPP')
MSs = lib_ms.AllMSs( glob.glob('mss-dd/TC*[0-9].MS'), s )
       
logger.info('Add columns...')
MSs.run('addcol2ms.py -m $pathMS -c CORRECTED_DATA,SUBTRACTED_DATA', log='$nameMS_addcol.log', commandType='python')

###############################################################
logger.info('BL-based smoothing...')
#MSs.run('BLsmooth.py -f 1.0 -r -i DATA -o SMOOTHED_DATA $pathMS', log='$nameMS_smooth.log', commandType='python')

# setup initial model
mosaic_image = lib_img.Image(sorted(glob.glob('self/images/wideM-[0-9]-MFS-image.fits'))[-1], userReg = userReg)
mosaic_image.selectCC()
rms_noise_pre = np.inf

for c in xrange(maxniter):
    logger.info('Starting cycle: %i' % c)

    lib_util.check_rm('img')
    os.makedirs('img')
    os.makedirs('ddcal/images/c%02i/regions' % c)

    ### group into patches of similar flux
    lsm = lsmtool.load(mosaic_image.skymodel_cut)
    lsm.group('tessellate', targetFlux='20Jy', root='Dir', applyBeam=False, method = 'wmean', pad_index=True)
    lsm.setPatchPositions(method='wmean')
    directions = lsm.getPatchPositions()
    patches = lsm.getPatchNames()
    logger.info("Created %i directions." % len(patches))

    # write file
    skymodel_cl = 'ddcal/skymodels/skymodel%02i_cluster.txt' % c
    lsm.write(skymodel_cl, format='makesourcedb', clobber=True)
    skymodel_cl_plot = 'ddcal/skymodels/skymodel%02i_cluster.png' % c
    lsm.plot(fileName=skymodel_cl_plot, labelBy='patch')

    # convert to blob
    skymodel_cl_skydb = skymodel_cl.replace('.txt','.skydb')
    lib_util.check_rm(skymodel_cl_skydb)
    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_cl, skymodel_cl_skydb), log='makesourcedb_cl.log', commandType='general' )
    s.run(check=True)

    ### create regions (using cluster directions)
    logger.info("Create regions.")
    lib_dd.make_voronoi_reg(directions, mosaic_image.maskname, outdir_reg='ddcal/images/c%02i/regions/' % c, out_mask='ddcal/masks/facets%02i.fits' % c, beam_reg='', png='ddcal/skymodels/voronoi%02i.png' % c)
    lsm.group('facet', facet='ddcal/masks/facets%02i.fits' % c)
    sizes = lsm.getPatchSizes(units='degree')

    # write file
    skymodel_voro = 'ddcal/skymodels/skymodel%02i_voro.txt' % c
    lsm.write(skymodel_voro, format='makesourcedb', clobber=True)
    skymodel_voro_plot = 'ddcal/skymodels/skymodel%02i_voro.png' % c
    lsm.plot(fileName=skymodel_voro_plot, labelBy='patch')

    # convert to blob
    skymodel_voro_skydb = skymodel_voro.replace('.txt','.skydb')
    lib_util.check_rm(skymodel_voro_skydb)
    s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_voro, skymodel_voro_skydb), log='makesourcedb_voro.log', commandType='general')
    s.run(check=True)

    del lsm

    ################################################################
    # Calibration
    logger.info('Calibrating...')
    MSs.run('DPPP '+parset_dir+'/DPPP-solDD.parset msin=$pathMS ddecal.h5parm=$pathMS/cal-c'+str(c)+'.h5 ddecal.sourcedb='+skymodel_cl_skydb, \
            log='$nameMS_solDD-c'+str(c)+'.log', commandType='DPPP')

    # Plot solutions
    logger.info('Running losoto...')
    lib_util.run_losoto(s, 'c'+str(c), [MS+'/cal-c'+str(c)+'.h5' for MS in MSs.getListStr()], [parset_dir+'/losoto-plot.parset'])
    os.system('mv plots-c'+str(c)+'* ddcal/plots')

    ############################################################
    # Empty the dataset
    logger.info('Set SUBTRACTED_DATA = DATA...')
    MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA"', log='$nameMS_taql1-c'+str(c)+'.log', commandType='general')

    logger.info('Subtraction...')
    MSs.run('DPPP '+parset_dir+'/DPPP-sub.parset msin=$pathMS sub.parmdb=$pathMS/cal-c'+str(c)+'.h5 sub.sourcedb='+skymodel_voro_skydb, \
                   log='$nameMS_sub-c'+str(c)+'.log', commandType='DPPP')

#    for i, p in enumerate(patches):
#        # predict - ms:MODEL_DATA
#        logger.info('Patch '+p+': predict...')
#        #pre.applycal.h5parm='+ms+'/cal-c'+str(c)+'.h5 pre.applycal.direction='+p, \
#        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+p, \
#                   log='$nameMS_pre1-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')
#
#        # corrupt - ms:MODEL_DATA -> ms:MODEL_DATA
#        logger.info('Patch '+p+': corrupt...')
#        MSs.run('DPPP '+parset_dir+'/DPPP-corrupt.parset msin=$pathMS cor1.parmdb=$pathMS/cal-c'+str(c)+'.h5 cor1.direction=['+p+'] \
#                                                                                 cor2.parmdb=$pathMS/cal-c'+str(c)+'.h5 cor2.direction=['+p+']', \
#                 log='$nameMS_corrupt1-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')
#
#        logger.info('Patch '+p+': subtract...')
#        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA"', log='$nameMS_taql2-c'+str(c)+'-p'+str(p)+'.log', commandType='general')

    ##############################################################
    # Imaging
    logger.info('Imaging...')

    ## TODO: test
    #logger.info('Empty imaging')
    #s.add('wsclean -reorder -name img/empty-c'+str(c)+' -datacolumn SUBTRACTED_DATA -size 3000 3000 \
    #        -mem 90 -j '+str(s.max_processors)+' -baseline-averaging 2.0 \
    #        -scale 10arcsec -weight briggs 0.0 -niter 0 -no-update-model-required -mgain 1 -pol I '+' '.join(mss), \
    #        log='wscleanEmpty-c'+str(c)+'.log', commandType='wsclean', processors='max')
    #s.run(check=True)

    for i, p in enumerate(patches):

        # add back single path - ms:SUBTRACTED_DATA -> ms:CORRECTED_DATA
        logger.info('Patch '+p+': add back...')
        MSs.run('DPPP '+parset_dir+'/DPPP-add.parset msin=$pathMS add.parmdb=$pathMS/cal-c'+str(c)+'.h5 add.sourcedb='+skymodel_voro_skydb+' add.directions=['+p+']', \
                   log='$nameMS_add-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')

#        #pre.applycal.h5parm='+ms+'/cal-c'+str(c)+'.h5 pre.applycal.direction='+p, \
#        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+p, \
#                   log='$nameMS_pre2-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')
#
#        # corrupt - ms:MODEL_DATA -> ms:MODEL_DATA
#        logger.info('Patch '+p+': corrupt...')
#        MSs.run('DPPP '+parset_dir+'/DPPP-corrupt.parset msin=$pathMS cor1.parmdb=$pathMS/cal-c'+str(c)+'.h5 cor1.direction=['+p+'] cor2.parmdb=$pathMS/cal-c'+str(c)+'.h5 cor2.direction=['+p+']', \
#                 log='$nameMS_corrupt2-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')
#
#        logger.info('Patch '+p+': add...')
#        MSs.run('taql "update $pathMS set CORRECTED_DATA = SUBTRACTED_DATA + MODEL_DATA"', log='$nameMS_taql2-c'+str(c)+'-p'+str(p)+'.log', commandType='general')

        ### TEST
        #logger.info('Patch '+p+': phase shift and avg...')
        #lib_util.check_rm('mss_dd')
        #os.makedirs('mss_dd')
        #for ms in mss:
        #    msout = 'mss_dd/'+os.path.basename(ms)
        #    phasecentre = directions_shifts[p]
        #    s.add('DPPP '+parset_dir+'/DPPP-shiftavg.parset msin='+ms+' msout='+msout+' shift.phasecenter=['+str(phasecentre[0].degree)+'deg,'+str(phasecentre[1].degree)+'deg\]', \
        #        log=ms+'_shift-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')
        #s.run(check=True)
        #logger.info('Patch '+p+': corrupted image...')
        #clean('uncor-'+p, glob.glob('mss_dd/*MS'), size=sizes[i])
        ### end TEST

        # DD-correct - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
        logger.info('Patch '+p+': correct...')
        MSs.run('DPPP '+parset_dir+'/DPPP-cor.parset msin=$pathMS cor1.parmdb=$pathMS/cal-c'+str(c)+'.h5 cor1.direction=['+p+']', \
               log='$nameMS_cor-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')

        logger.info('Patch '+p+': phase shift and avg...')
        lib_util.check_rm('mss-dir')
        os.makedirs('mss-dir')
        phasecentre = directions[p]
        MSs.run('DPPP '+parset_dir+'/DPPP-avg.parset msin=$pathMS msout=mss-dir/$nameMS.MS msin.datacolumn=CORRECTED_DATA \
                shift.phasecenter=['+str(phasecentre[0].degree)+'deg,'+str(phasecentre[1].degree)+'deg\]', \
                log='$nameMS_shift-c'+str(c)+'-p'+str(p)+'.log', commandType='DPPP')
        
        logger.info('Patch '+p+': imaging...')
        clean(p, lib_ms.AllMSs('mss-dir/*MS'), size=sizes[i])

    ##############################################################
    # Mosaiching

    directions = []
    for image, region in zip( sorted(glob.glob('img/ddcalM-Dir*MFS-image.fits')), sorted(glob.glob('ddcal/images/c%02i/regions/Dir*' % c)) ):
        directions.append( Image(image, regionFacet = region, user_mask = user_mask) )

    logger.info('Mosaic: image...')
    images = ' '.join([image.imagename for image in directions])
    masks = ' '.join([image.regionFacet for image in directions])
    mosaic_imagename = 'img/mos-MFS-image.fits'
    s.add('mosaic.py --image '+images+' --masks '+masks+' --output '+mosaic_imagename, log='mosaic-img-c'+str(c)+'.log', commandType='python')
    s.run(check=True)

    logger.info('Mosaic: residuals...')
    images = ' '.join([image.imagename.replace('image', 'residual') for image in directions])
    masks = ' '.join([image.regionFacet for image in directions])
    mosaic_residual = 'img/mos-MFS-residual.fits'
    s.add('mosaic.py --image '+images+' --masks '+masks+' --output '+mosaic_residual, log='mosaic-res-c'+str(c)+'.log', commandType='python')
    s.run(check=True)

    # prepare new skymodel
    skymodels = []
    for image in directions:
        image.select_cc()
        skymodels.append(image.skymodel_cut)
    lsm = lsmtool.load(skymodels[0])
    for skymodel in skymodels[1:]:
        lsm2 = lsmtool.load(skymodel)
        lsm.concatenate(lsm2)
    lsm.write('ddcal/images/c%02i/mos-sources-cut.txt' % c, format='makesourcedb', clobber=True)

    os.system('cp img/*M*MFS-image.fits img/mos-MFS-image.fits img/mos-MFS-residual.fits ddcal/images/c%02i' % c )
    mosaic_image = Image('ddcal/images/c%02i/mos-MFS-image.fits' % c, user_mask = user_mask)

    # get noise, if larger than 95% of prev cycle: break
    rms_noise = lib_img.Image(mosaic_residual).getNoise()
    logger.info('RMS noise: %f' % rms_noise)
    if rms_noise > 0.95 * rms_noise_pre: break
    rms_noise_pre = rms_noise
