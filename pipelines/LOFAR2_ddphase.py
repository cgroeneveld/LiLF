#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Pipeline for direction dependent calibration

import sys, os, glob, re, pickle
import numpy as np
import pyrap.tables as pt
import lsmtool

#######################################################
from LiLF import lib_ms, lib_img, lib_util, lib_log, lib_dd_parallel
logger_obj = lib_log.Logger('pipeline-dd.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir = logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-dd.walker')

# parse parset
parset = lib_util.getParset()
parset_dir = parset.get('LOFAR_dd-parallel','parset_dir')
maxniter = parset.getint('LOFAR_dd-parallel','maxniter')
calFlux = parset.getfloat('LOFAR_dd-parallel','calFlux')
userReg = parset.get('model','userReg')
aterm_imaging = False

MSs_self = lib_ms.AllMSs( glob.glob('mss/TC*[0-9].MS'), s )


# make beam
phasecentre = MSs_self[0].getPhaseCentre()
fwhm = MSs_self[0].getFWHM(freq='mid')
# for frequency scaling
freqscale = np.mean(MSs_self.getFreqs())/58.e6
smoothfactor = 0.002 if MSs_self.isHBA else 0.01

def clean(p, MSs, size, res='normal', apply_beam=False):
    """
    p = patch name
    mss = list of mss to clean
    size = in deg of the image
    """

    # set pixscale and imsize
    pixscale = MSs[0].getResolution()
    if res == 'normal':
        pixscale = float('%.1f'%(pixscale/2.5))
    elif res == 'high':
        pixscale = float('%.1f'%(pixscale/3.5))
    elif res == 'low':
        pass # no change

    imsize = [0,0]
    imsize[0] = int(size[0]*1.1/(pixscale/3600.)) # add 10%
    imsize[1] = int(size[1]*1.1/(pixscale/3600.)) # add 10%
    imsize[0] += imsize[0] % 2
    imsize[1] += imsize[1] % 2
    if imsize[0] < 64: imsize[0] == 64
    if imsize[1] < 64: imsize[1] == 64

    logger.info('Image size: '+str(imsize)+' - Pixel scale: '+str(pixscale))

    maxuv_l = None
    if res == 'normal':
        weight = 'briggs -0.1'
    elif res == 'high':
        weight = 'briggs -0.6'
    elif res == 'low':
        weight = 'briggs 0'
        maxuv_l = 3500*freqscale
    channels_out = 9 if MSs.isLBA else 6
    # CLEAN
    logger.info('Cleaning ('+str(p)+')...')
    imagename = 'img/ddcal-'+str(p)
    if apply_beam:
        # idg cannot be used for rectangular images.
        # multiscale:
        imsize = [max(imsize),max(imsize)]
        lib_util.run_wsclean(s, 'wsclean-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, size=imsize, save_source_list='', scale=str(pixscale)+'arcsec', \
                             weight=weight, niter=200000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
                             use_idg='', grid_with_beam='', use_differential_lofar_beam='', \
                             multiscale='', multiscale_scales='0,10,20', \
                             local_rms='', parallel_deconvolution=2048, auto_threshold=2.0, auto_mask=3.5, \
                             join_channels='', fit_spectral_pol=3, channels_out=channels_out, deconvolution_channels=3)
    else:
        lib_util.run_wsclean(s, 'wsclean-'+str(p)+'.log', MSs.getStrWsclean(), name=imagename, size=imsize, save_source_list='', scale=str(pixscale) + 'arcsec', \
                             weight=weight, niter=200000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, \
                             baseline_averaging='', parallel_gridding=2, parallel_deconvolution=2048, auto_threshold=1.0, \
                             local_rms='', nmiter=4, multiscale='', multiscale_scales='0,7,16,27,40',\
                             auto_mask=3.0, join_channels='', fit_spectral_pol=3, channels_out=channels_out, deconvolution_channels=3)

    os.system('cat logs/wsclean-'+str(p)+'.log| grep "background noise"')


#############################################################
with w.if_todo('cleaning'):
    logger.info('Cleaning...')
    lib_util.check_rm('ddcal')
    os.makedirs('ddcal/masks')
    os.makedirs('ddcal/plots')
    os.makedirs('ddcal/images')
    os.makedirs('ddcal/solutions')
    os.makedirs('ddcal/skymodels')
### DONE

############################################################
# use SUBTRACTED_DATA (no pre-correction - subtraction would not work) or CORRECTED_DATA (DIE iono correction)?
logger.info('Copy data...')
if not os.path.exists('mss-dd'):
    os.makedirs('mss-dd')
    MSs_self.run('DPPP '+parset_dir+'/DPPP-avg.parset msin=$pathMS msout=mss-dd/$nameMS.MS msin.datacolumn=CORRECTED_DATA avg.freqstep=1 avg.timestep=1', \
                log='$nameMS_avg.log', commandType='DPPP')
MSs = lib_ms.AllMSs( glob.glob('mss-dd/TC*[0-9].MS'), s )
       
with w.if_todo('addcol'):
    logger.info('Add columns...')
    MSs.run('addcol2ms.py -m $pathMS -c CORRECTED_DATA,SUBTRACTED_DATA -i DATA', log='$nameMS_addcol.log', commandType='python')

##############################################################
# setup initial model
MSs[0].makeBeamReg('ddcal/beam.reg', freq='mid')

if MSs.isLBA: # For LOFAR2.0 sim use larger beam for HBA
    beamReg = 'ddcal/beam.reg'
elif MSs.isHBA:
    beamReg = 'ddcal/lbabeam.reg'
mosaic_image = lib_img.Image(sorted(glob.glob('self/images/wideM-[0-9]-MFS-image.fits'))[-1], beamReg=beamReg, userReg = userReg)
logger.warning('EDITED for HBA phases')
# mosaic_image.makeMask() # create mask of extended bright stuff wideM-0-mask.fits
# mosaic_image.selectCC() # create wideM-o-sources-cut.skydb
# lsm = lsmtool.load(mosaic_image.skymodel_cut)

rms_noise_pre = np.inf
dir_from_lba = True

for c in range(maxniter):
    logger.info('Starting cycle: %i' % c)
    if c>=1: directions_old = directions

    with w.if_todo('delimg-c%02i' % c):
        lib_util.check_rm('img')
        os.makedirs('img')
    ### DONE

    skymodel_cl = 'ddcal/skymodels/skymodel%02i_cluster.txt' % c
    skymodel_cl_skydb = skymodel_cl.replace('.txt','.skydb')
    skymodel_rest = 'ddcal/skymodels/skymodel%02i_rest.txt' % c
    skymodel_rest_skydb = skymodel_rest.replace('.txt','.skydb')
    skymodel_voro = 'ddcal/skymodels/skymodel%02i_voro.txt' % c
    skymodel_voro_skydb = skymodel_voro.replace('.txt','.skydb')

    picklefile = 'ddcal/directions-c%02i.pickle' % c
    if dir_from_lba:
        directions = []

        if not os.path.exists('ddcal/masks/regions-c%02i' % c): os.makedirs('ddcal/masks/regions-c%02i' % c)
        if not os.path.exists('ddcal/images/c%02i' % c): os.makedirs('ddcal/images/c%02i' % c)
        mask_voro = 'ddcal/masks/facets%02i.fits' % c

        lsm = lsmtool.load(skymodel_cl)
        for name, flux in zip(lsm.getPatchNames(), lsm.getColValues('I', aggregate='sum')):
            direction = lib_dd_parallel.Direction(name)
            position = [ lsm.getPatchPositions()[name][0].deg, lsm.getPatchPositions()[name][1].deg ]
            direction.set_position( position, cal=True )
            direction.set_flux(flux, cal=True)
            directions.append(direction)

        lib_dd_parallel.make_voronoi_reg(directions, mosaic_image.maskname, \
                                         outdir_reg='ddcal/masks/regions-c%02i' % c, out_mask=mask_voro,
                                         png='ddcal/masks/voronoi%02i.png' % c)
        [d.add_mask_voro(mask_voro) for d in directions]
        pickle.dump(directions, open(picklefile, "wb"))

    elif not os.path.exists(picklefile):
        directions = []

        if not os.path.exists('ddcal/masks/regions-c%02i' % c): os.makedirs('ddcal/masks/regions-c%02i' % c)
        if not os.path.exists('ddcal/images/c%02i' % c): os.makedirs('ddcal/images/c%02i' % c)
        mask_voro = 'ddcal/masks/facets%02i.fits' % c

        ### group into patches corresponding to the mask islands
        mask_cl = mosaic_image.imagename.replace('image.fits', 'mask-cl.fits')
        # this mask is with no user region, done to isolate only bight compact sources
        if not os.path.exists(mask_cl): 
            mosaic_image.beamReg = 'ddcal/beam.reg'
            if MSs.isLBA:
                if c == 0:
                    mosaic_image.makeMask(threshisl=7, atrous_do=False, remove_extended_cutoff=0.001, maskname=mask_cl, only_beam=True)
                else:
                    mosaic_image.makeMask(threshisl=5, atrous_do=False, remove_extended_cutoff=0.01, maskname=mask_cl, only_beam=True)
            elif MSs.isHBA:
                mosaic_image.makeMask(threshisl=7, atrous_do=False, remove_extended_cutoff=0.0001, maskname=mask_cl, only_beam=True)
            else:
                logger.error('No LBA or HBA.')

        lsm = lsmtool.load(mosaic_image.skymodel_cut)

        lsm.group(mask_cl, root='Isl')
        # this removes all sources not in the mask-cl
        lsm.select('Patch = Isl.*', useRegEx=True)
        lsm.info()
        # this regroup sources
        x = lsm.getColValues('RA',aggregate='wmean')
        y = lsm.getColValues('Dec',aggregate='wmean')
        flux = lsm.getColValues('I',aggregate='sum')
        if MSs.isLBA:
            grouper = lib_dd_parallel.Grouper(list(zip(x,y)),flux,look_distance=0.4,kernel_size=0.3,grouping_distance=0.1)
        elif MSs.isHBA:
            grouper = lib_dd_parallel.Grouper(list(zip(x, y)), flux, look_distance=0.3, kernel_size=0.2, grouping_distance=0.1)
        grouper.run()
        clusters = grouper.grouping()
        grouper.plot()
        os.system('mv grouping*png ddcal/plots/')
        patchNames = lsm.getPatchNames()
    
        logger.info('Merging nearby sources...')
        for cluster in clusters:
            patches = patchNames[cluster]
            print('merging:', cluster, patches)
            if len(patches) > 1:
                lsm.merge(patches.tolist())
        print(lsm.info())
        print(lsm.getColValues('I',aggregate='sum'))
        lsm.select('I >= %f Jy' % float(calFlux), aggregate='sum')
        print(lsm)
    
        # keep track of CC names used for calibrators so not to subtract them afterwards
        cal_names = lsm.getColValues('Name')
    
        lsm.setPatchPositions(method='wmean') # calculate patch weighted centre for tassellation
        for name, flux in zip(lsm.getPatchNames(), lsm.getColValues('I', aggregate='sum')):
            direction = lib_dd_parallel.Direction(name)
            position = [ lsm.getPatchPositions()[name][0].deg, lsm.getPatchPositions()[name][1].deg ]
            direction.set_position( position, cal=True )
            direction.set_flux(flux, cal=True)
            directions.append(direction)


        # write file
        lsm.write(skymodel_cl, format='makesourcedb', clobber=True)
        skymodel_cl_plot = 'ddcal/masks/skymodel%02i_cluster.png' % c
        lsm.plot(fileName=skymodel_cl_plot, labelBy='patch')
        lsm.setColValues('name', [x.split('_')[-1] for x in lsm.getColValues('patch')]) # just for the region - this makes this lsm useless
        lsm.write('ddcal/masks/regions-c%02i/cluster.reg' % c, format='ds9', clobber=True)
        del lsm
    
        # convert to blob
        lib_util.check_rm(skymodel_cl_skydb)
        s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_cl, skymodel_cl_skydb), log='makesourcedb_cl.log', commandType='general' )
        s.run(check=True)
        
        ### select the rest of the sources to be subtracted
        lsm = lsmtool.load(mosaic_image.skymodel_cut)
        names = lsm.getColValues('Name')
        lsm.remove( np.array([ i for i, name in enumerate(names) if name in cal_names ]) )
        lsm.ungroup()
        logger.info("Total flux in rest field %i Jy" % np.sum(lsm.getColValues('I')) )
        
        # write file
        lsm.write(skymodel_rest, format='makesourcedb', clobber=True)
        skymodel_rest_plot = 'ddcal/masks/skymodel%02i_rest.png' % c
        lsm.plot(fileName=skymodel_rest_plot, labelBy='patch')
           
        # convert to blob
        lib_util.check_rm(skymodel_rest_skydb)
        s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_rest, skymodel_rest_skydb), log='makesourcedb_rest.log', commandType='general')
        s.run(check=True)
        
        ### create regions (using cluster directions)
        logger.info("Create regions.")
        lsm = lsmtool.load(mosaic_image.skymodel_cut)

        lib_dd_parallel.make_voronoi_reg(directions, mosaic_image.maskname, \
                outdir_reg='ddcal/masks/regions-c%02i' % c, out_mask=mask_voro, png='ddcal/masks/voronoi%02i.png' % c)
        lsm.group('facet', facet=mask_voro, root='Isl_patch')
        [ d.add_mask_voro(mask_voro) for d in directions ]
    
        # write file
        lsm.write(skymodel_voro, format='makesourcedb', clobber=True)
        skymodel_voro_plot = 'ddcal/masks/skymodel%02i_voro.png' % c
        lsm.plot(fileName=skymodel_voro_plot, labelBy='patch')
        lsm.setColValues('name', [x.split('_')[-1] for x in lsm.getColValues('patch')]) # just for the region - this makes this lsm useless
        lsm.write('ddcal/masks/regions-c%02i/voro.reg' % c, format='ds9', clobber=True)
        del lsm
    
        # convert to blob
        lib_util.check_rm(skymodel_voro_skydb)
        s.add('makesourcedb outtype="blob" format="<" in="%s" out="%s"' % (skymodel_voro, skymodel_voro_skydb), log='makesourcedb_voro.log', commandType='general')
        s.run(check=True)

        pickle.dump( directions, open( picklefile, "wb" ) )

    else:
        directions = pickle.load( open( picklefile, "rb" ) )
        mask_voro = 'ddcal/masks/facets%02i.fits' % c

    logger.info("Created %i bright sources" % len(directions))
    tot_flux = np.sum([d.flux_cal for d in directions])
    logger.info("Total flux of bright sources %i Jy" % tot_flux)
    
    logger.debug("Islands' info:")
    for i, d in enumerate(directions):
        logger.info("%s: Flux=%f (coord: %s - size: %s deg)" % ( d.name, d.flux_cal, str(d.position_cal), str(d.size_facet) ) )
    ###############################################################
    with w.if_todo('subtract-faint-c0{}'.format(c)):
        logger.info('Subtraction rest_field...')
        # Predict - ms:MODEL_DATA
        logger.info('Add rest_field to MODEL_DATA...')
        MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_rest_skydb,log='$nameMS_pre-c'+str(c)+'.log', commandType='DPPP')
        # Empty dataset from faint sources
        logger.info('Set SUBTRACTED_DATA = DATA - MODEL_DATA...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
        # Smoothing - ms:SUBTRACTED_DATA -> ms:SMOOTHED_DATA
        logger.info('BL-based smoothing...')
        MSs.run('BLsmooth.py -c 8 -f {} -r -i SUBTRACTED_DATA -o SMOOTHED_DATA $pathMS'.format(smoothfactor),
                log='$nameMS_smooth-c'+str(c)+'.log', commandType='python')

    # Calibrate
    with w.if_todo('calibrate-c%02i' % c):
        # Calibration - ms:SMOOTHED_DATA
        logger.info('Phase calibration...')
        dirs = '[Isl_patch_12,Isl_patch_31,Isl_patch_10,Isl_patch_35]'
        if MSs.isLBA:
            MSs.run('DPPP ' + parset_dir + '/DPPP-solPh.parset msin=$pathMS \
                     sol.h5parm=$pathMS/cal-ph-c' + str(c) + '.h5 sol.sourcedb=' + skymodel_cl_skydb \
                    +' sol.antennaconstraint=[[CS001LBA,CS002LBA,CS003LBA,CS004LBA,CS005LBA,CS006LBA,'
                     'CS007LBA,CS011LBA,CS013LBA,CS017LBA,CS021LBA,CS024LBA,CS026LBA,CS028LBA,CS030LBA,'
                     'CS031LBA,CS032LBA,CS101LBA,CS103LBA,CS201LBA,CS301LBA,CS302LBA,CS401LBA,CS501LBA]] '
                     'sol.directions={} sol.uvlambdamin={} '.format(dirs,100*freqscale),
                    log='$nameMS_solPh-c' + str(c) + '.log', commandType='DPPP')

            lib_util.run_losoto(s, 'ph-c'+str(c), [ms+'/cal-ph-c'+str(c)+'.h5' for ms in MSs.getListStr()], \
                                [parset_dir+'/losoto-plot-ph.parset'])
        elif MSs.isHBA:
            MSs.run('DPPP '+parset_dir+'/DPPP-solPh.parset msin=$pathMS \
                    sol.h5parm=$pathMS/cal-ph-c'+str(c)+'.h5 sol.sourcedb='+skymodel_cl_skydb
                    +' sol.antennaconstraint=[[CS001HBA0,CS001HBA1,'
                                             'CS002HBA0,CS002HBA1,'
                                             'CS003HBA0,CS003HBA1,'
                                             'CS004HBA0,CS004HBA1,'
                                             'CS005HBA0,CS005HBA1,'
                                             'CS006HBA0,CS006HBA1,'
                                             'CS007HBA0,CS007HBA1,'
                                             'CS011HBA0,CS011HBA1,'
                                             'CS013HBA0,CS013HBA1,'
                                             'CS017HBA0,CS017HBA1,'
                                             'CS021HBA0,CS021HBA1,'
                                             'CS024HBA0,CS024HBA1,'
                                             'CS026HBA0,CS026HBA1,'
                                             'CS028HBA0,CS028HBA1,'
                                             'CS030HBA0,CS030HBA1,'
                                             'CS031HBA0,CS031HBA1,'
                                             'CS032HBA0,CS032HBA1,'
                                             'CS101HBA0,CS101HBA1,'
                                             'CS103HBA0,CS103HBA1,'
                                             'CS201HBA0,CS201HBA1,'
                                             'CS301HBA0,CS301HBA1,'
                                             'CS302HBA0,CS302HBA1,'
                                             'CS401HBA0,CS401HBA1,'
                                             'CS501HBA0,CS501HBA1]] '
                     'sol.directions={} sol.uvlambdamin={} '.format(dirs, 100 * freqscale),
                    log='$nameMS_solPh-c'+str(c)+'.log', commandType='DPPP')
            lib_util.run_losoto(s, 'ph-c' + str(c), [ms + '/cal-ph-c' + str(c) + '.h5' for ms in MSs.getListStr()], \
                                [parset_dir + '/losoto-plot-ph.parset'])

        # Plot solutions
        os.system('mv plots-ph-c' + str(c) + ' ddcal/plots')
        os.system('mv cal-ph-c' + str(c) + '.h5 ddcal/solutions')
    ### DONE
    sys.exit()

    ###########################################################
    # use idg and A-term to correct the data, single image
    if aterm_imaging:

        #wsclean -mem 90.0 -scale 0.0004166666666666667 -aterm-config /beegfs/rafferty/Data/LOFAR/Screens/Factor_sim/pipelines/image_1/sector_3/chunk9.ms.make_aterm_config -multiscale-scales 0 -size 1500 1500 -deconvolution-channels 4 -fits-mask /beegfs/rafferty/Data/LOFAR/Screens/Factor_sim/pipelines/image_1/sector_3/chunk9.ms.premask -j 6 -auto-mask 3.6 -idg-mode hybrid -channels-out 12 -local-rms-window 50 -mgain 0.5 -minuv-l 80.0 -fit-spectral-pol 3 -maxuv-l 1000000.0 -weighting-rank-filter 3 -aterm-kernel-size 32 -temp-dir /tmp -name /beegfs/rafferty/Data/LOFAR/Screens/Factor_sim/pipelines/image_1/sector_3/chunk9.ms.image -padding 1.2 -pol I -multiscale-shape gaussian -auto-threshold 1.0 -local-rms-method rms-with-min -weight briggs -0.5 -niter 13635 -no-update-model-required -multiscale -fit-beam -reorder -save-source-list -local-rms -join-channels -use-idg -apply-primary-beam -nmiter 4
        lib_util.run_wsclean(s, 'wscleanDD-c%02i.log' %c, MSs.getStrWsclean(), name='img/wideDD-c%02i' %c, save_source_list='', size=imsize, scale=str(pixscale)+'arcsec', \
            weight=weight, niter=100000, no_update_model_required='', minuv_l=30, maxuv_l=maxuv_l, mgain=0.85, \
            use_idg='', grid_with_beam='', use_differential_lofar_beam='', beam_aterm_update=400, \
            multiscale='', multiscale_scales='0,10,20,40,80', \
            parallel_deconvolution=512, local_rms='', auto_threshold=0.75, auto_mask=1.5, fits_mask=im.maskname, \
            join_channels='', fit_spectral_pol=3, channels_out=9, deconvolution_channels=3)

        # TODO: put proper names
        os.system('cp img/wideDD-c%02i.MFS-image.fits ddcal/images/c%02i' % (c,c) )
        mosaic_image = lib_img.Image('ddcal/images/c%02i/mos-MFS-image.fits' % c, userReg = userReg)

    ###########################################################
    # facet imaging
    ###########################################################
    # Subtraction
    with w.if_todo('empty-c%02i' % c):
        logger.info('Subtraction...')
        # Copy DATA -> SUBTRACTED_DATA
        logger.info('Set SUBTRACTED_DATA = DATA...')
        MSs.run('taql "update $pathMS set SUBTRACTED_DATA = DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
        for i, d in enumerate(directions):
            # predict - ms:MODEL_DATA
            logger.info('Patch '+d.name+': predict...')
            MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+d.name, \
                    log='$nameMS_pre1-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

            # corrupt TEC - ms:MODEL_DATA -> ms:MODEL_DATA
            logger.info('Patch '+d.name+': corrupt...')
            MSs.run('DPPP '+parset_dir+'/DPPP-corrupt1.parset msin=$pathMS \
                    cor.parmdb=ddcal/solutions/cal-tec-c'+str(c)+'.h5 cor.correction=tec000 cor.direction=['+d.name+']', \
                    log='$nameMS_corrupt1-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')
            # if c>0:
            #     MSs.run('DPPP '+parset_dir+'/DPPP-corrupt1.parset msin=$pathMS \
            #         cor.parmdb=ddcal/solutions/cal-g-c'+str(c)+'.h5 cor.correction=amplitude000 cor.direction=['+d.name+']', \
            #         log='$nameMS_corrupt1-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

            logger.info('Patch '+d.name+': subtract...')
            MSs.run('taql "update $pathMS set SUBTRACTED_DATA = SUBTRACTED_DATA - MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'-'+d.name+'.log', commandType='general')
    ### DONE

    ### TESTTESTTEST: empty image
    #MSs.run('taql "update $pathMS set CORRECTED_DATA = SUBTRACTED_DATA"', log='$nameMS_taql-c'+str(c)+'.log', commandType='general')
    # clean('empty-c'+str(c), MSs, size=(fwhm,fwhm), res='low', data_column='SUBTRACTED_DATA')
    ###

    ###########################################################
    # Add back
    logger.info('Facet imaging...')
    for i, d in enumerate(directions):

        with w.if_todo('facet-%s-c%02i' % (d.name,c)):
            #TODO: see if we can phase shift and average before predict-corrupt=correct
            # predict - ms:MODEL_DATA
            logger.info('Patch '+d.name+': predict...')
            MSs.run('DPPP '+parset_dir+'/DPPP-predict.parset msin=$pathMS pre.sourcedb='+skymodel_voro_skydb+' pre.sources='+d.name, \
                       log='$nameMS_pre2-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

            # # corrupt G - ms:MODEL_DATA -> ms:MODEL_DATA
            # logger.info('Patch '+d.name+': corrupt...')
            # MSs.run('DPPP '+parset_dir+'/DPPP-corrupt1.parset msin=$pathMS \
            #         cor.parmdb=ddcal/solutions/cal-g-c'+str(c)+'.h5 cor.correction=phase000 cor.direction=['+d.name+']', \
            #         log='$nameMS_corrupt2-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')
            # corrupt tec - ms:MODEL_DATA -> ms:MODEL_DATA
            logger.info('Patch '+d.name+': corrupt...')
            MSs.run('DPPP '+parset_dir+'/DPPP-corrupt1.parset msin=$pathMS \
                    cor.parmdb=ddcal/solutions/cal-tec-c'+str(c)+'.h5 cor.correction=tec000 cor.direction=['+d.name+']', \
                    log='$nameMS_corrupt2-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

            logger.info('Patch '+d.name+': add...')
            MSs.run('taql "update $pathMS set CORRECTED_DATA = SUBTRACTED_DATA + MODEL_DATA"', log='$nameMS_taql-c'+str(c)+'-'+d.name+'.log', commandType='general')

            # correct G - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
            # logger.info('Patch '+d.name+': correct...')
            # MSs.run('DPPP '+parset_dir+'/DPPP-correct1.parset msin=$pathMS \
            #         cor.parmdb=ddcal/solutions/cal-g-c'+str(c)+'.h5 cor.correction=phase000 cor.direction=['+d.name+']', \
            #         log='$nameMS_correct-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')
            # correct TEC - ms:CORRECTED_DATA -> ms:CORRECTED_DATA
            logger.info('Patch '+d.name+': correct...')
            MSs.run('DPPP '+parset_dir+'/DPPP-correct1.parset msin=$pathMS \
                    cor.parmdb=ddcal/solutions/cal-tec-c'+str(c)+'.h5 cor.correction=tec000 cor.direction=['+d.name+']', \
                    log='$nameMS_correct-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

            # phase shift
            logger.info('Patch '+d.name+': phase shift and avg...')
            lib_util.check_rm('mss-dir')
            os.makedirs('mss-dir')
            MSs.run('DPPP '+parset_dir+'/DPPP-shiftavg.parset msin=$pathMS msout=mss-dir/$nameMS.MS msin.datacolumn=CORRECTED_DATA \
                    shift.phasecenter=['+str(d.position_facet[0])+'deg,'+str(d.position_facet[1])+'deg\]', \
                    log='$nameMS_shift-c'+str(c)+'-'+d.name+'.log', commandType='DPPP')

            logger.info('Patch '+d.name+': imaging...')
            clean(d.name, lib_ms.AllMSs( glob.glob('mss-dir/*MS'), s ), size=d.size_facet, apply_beam=False)
            # apply_beam: IDG does not yet support square images, so do not use idg...
        ### DONE

    ##############################################################
    # Mosaiching

    with w.if_todo('mosaic-c%02i' % c):
        # reorder in increasing isl_num order
        isl_nums = [d.isl_num for d in directions]
        directions = [d for _, d in sorted(zip(isl_nums,directions))]

        for d in directions:
            # Use pb corrected images
            d.image = lib_img.Image('img/ddcal-%s-MFS-image.fits' % d.name, userReg = userReg)
            d.image_res = lib_img.Image('img/ddcal-%s-MFS-residual.fits' % d.name, userReg = userReg)
            # d.image_low = lib_img.Image('img/ddcal-%s-low-MFS-image.fits' % d.name, userReg = userReg)
            # d.image_high = lib_img.Image('img/ddcal-%s-high-MFS-image.fits' % d.name, userReg = userReg)

            # restrict skymodel to facet
            d.image.makeMask(threshisl=5)
            d.image.selectCC(checkBeam=False) # checkBeam edited since lib_img was changed
            # try:
            try:
                lsm = lsmtool.load(d.image.skymodel_cut)
                lsm.group('facet', facet=mask_voro, root='Isl_patch' )
                lsm.select('Patch = Isl_patch_%i' % d.isl_num )
                lsm.write(d.image.skymodel_cut, format='makesourcedb', clobber=True)
            except:
                print(lsm.info())
                logger.error("No sources in facet %s?" % d.name)
                pass

        logger.info('Mosaic: image...')
        image_files = ' '.join([d.image.imagename for d in directions])
        mosaic_image_file = 'img/mos-MFS-image.fits'
        s.add('mosaic.py --image '+image_files+' --mask '+mask_voro+' --output '+mosaic_image_file, log='mosaic-img-c'+str(c)+'.log', commandType='python')
        s.run(check=True)

        logger.info('Mosaic: residual image...')
        image_files = ' '.join([d.image_res.imagename for d in directions])
        mosaic_residual_file = 'img/mos-MFS-residual.fits'
        s.add('mosaic.py --image '+image_files+' --mask '+mask_voro+' --output '+mosaic_residual_file, log='mosaic-res-c'+str(c)+'.log', commandType='python')
        s.run(check=True)

        # if c>=2:
        #     logger.info('Mosaic: low-res image...')
        #     image_files = ' '.join([d.image_low.imagename for d in directions])
        #     mosaic_image_low_file = 'img/mos-low-MFS-image.fits'
        #     s.add('mosaic.py --image '+image_files+' --mask '+mask_voro+' --output '+mosaic_image_low_file, log='mosaic-img-low-c'+str(c)+'.log', commandType='python')
        #     s.run(check=True)
        #
        #     logger.info('Mosaic: high-res image...')
        #     image_files = ' '.join([d.image_high.imagename for d in directions])
        #     mosaic_image_high_file = 'img/mos-high-MFS-image.fits'
        #     s.add('mosaic.py --image '+image_files+' --mask '+mask_voro+' --output '+mosaic_image_high_file, log='mosaic-img-high-c'+str(c)+'.log', commandType='python')
        #     s.run(check=True)

        # prepare new skymodel
        lsm = lsmtool.load(directions[0].image.skymodel_cut)
        lsm.ungroup()
        for image in [d.image for d in directions[1:]]:
            if os.path.exists(image.skymodel_cut):
                lsm2 = lsmtool.load(image.skymodel_cut)
                lsm2.ungroup()
                lsm.concatenate(lsm2, keep='all')
        lsm.write('ddcal/images/c%02i/mos-sources-cut.txt' % c, format='makesourcedb', clobber=True)
        del lsm

        # os.system('cp img/*M*MFS-image.fits img/mos*.fits ddcal/images/c%02i' % c )
        os.system('cp img/*M*MFS-image.fits img/mos*.fits ddcal/images/c%02i' % c )
    ### DONE

    mosaic_image = lib_img.Image('ddcal/images/c%02i/mos-MFS-image.fits' % c, userReg = userReg)
    #mosaic_image_pb = lib_img.Image('ddcal/images/c%02i/mos-MFS-image-pb.fits' % c, userReg = userReg)

    mosaic_image.makeMask(threshisl=3, atrous_do=True) # used in the faceting function
    # get noise, if larger than 95% of prev cycle: break
    rms_noise = mosaic_image.getNoise()
    logger.info('RMS noise: %f' % rms_noise)
    if rms_noise > rms_noise_pre: break
    rms_noise_pre = rms_noise
