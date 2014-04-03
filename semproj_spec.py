import os
import subprocess
from ciao_contrib.runtool import *


# Convert the event FITS file into a fluximage, in order to use it to plot the
# x-ray image
def img(ID):
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Data/'
             + ID + '/repro/')
    os.system('fluximage evt2_df_ccd7.fits aplpy_' + ID +
              ' bands=0.5:8:3 binsize=1')


# Create the spectrum for the source region
def spec(ID, double):
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Data/'
             + ID + '/repro/')
    if ID in double:
        os.system('punlearn specextract')
        os.system('pset specextract infile="evt2_df_ccd7.fits"'
                  + '[sky=region(source_1.reg)]"')
        os.system('pset specextract outroot='+ID+'_1_spec')
        os.system('pset specextract bkgfile="evt2_df_ccd7.fits'
                  + '[sky=region(back.reg)]"')
        os.system('pset specextract weight=no correctpsf=yes')
        os.system('specextract')
        os.system('punlearn specextract')
        os.system('pset specextract infile="evt2_df_ccd7.fits'
                  + '[sky=region(source_2.reg)]"')
        os.system('pset specextract outroot='+ID+'_2_spec')
        os.system('pset specextract bkgfile="evt2_df_ccd7.fits'
                  + '[sky=region(back.reg)]"')
        os.system('pset specextract weight=no correctpsf=yes')
        os.system('specextract')
    else:
        os.system('punlearn specextract')
        os.system('pset specextract infile="evt2_df_ccd7.fits'
                  + '[sky=region(source.reg)]"')
        os.system('pset specextract outroot='+ID+'_spec')
        os.system('pset specextract bkgfile="evt2_df_ccd7.fits'
                  + '[sky=region(back.reg)]"')
        os.system('pset specextract weight=no correctpsf=yes')
        os.system('specextract')


# Create the spectrum for the extended source region
def spec_ext(ID, double):
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Data/'
             + ID + '/repro/')
    if ID in double:
        os.system('punlearn specextract')
        os.system('pset specextract infile="evt2_df_ccd7.fits'
                  + '[sky=region(source_1_ext.reg)]"')
        os.system('pset specextract outroot='+ID+'_1_spec_ext')
        os.system('pset specextract bkgfile="evt2_df_ccd7.fits'
                  + '[sky=region(back.reg)]"')
        os.system('pset specextract weight=no correctpsf=yes')
        os.system('specextract')
        os.system('punlearn specextract')
        os.system('pset specextract infile="evt2_df_ccd7.fits'
                  + '[sky=region(source_2_ext.reg)]"')
        os.system('pset specextract outroot='+ID+'_2_spec_ext')
        os.system('pset specextract bkgfile="evt2_df_ccd7.fits'
                  + '[sky=region(back.reg)]"')
        os.system('pset specextract weight=no correctpsf=yes')
        os.system('specextract')
    else:
        os.system('punlearn specextract')
        os.system('pset specextract infile="evt2_df_ccd7.fits'
                  + '[sky=region(source_ext.reg)]"')
        os.system('pset specextract outroot='+ID+'_spec_ext')
        os.system('pset specextract bkgfile="evt2_df_ccd7.fits'
                  + '[sky=region(back.reg)]"')
        os.system('pset specextract weight=no correctpsf=yes')
        os.system('specextract')


# Run the deflare script on the event file for a single source, excluding the
# source region
def deflareobs(evt_file, fov_file):
    dmcopy(infile=evt_file + '[energy=500:7000,ccd_id=7,'
           + 'sky=region('+fov_file+'[ccd_id=7])]', outfile='evt2_c7.fits',
           clobber='yes')
    print 'hit return to exit deflare plot'
    dmextract(infile='evt2_c7.fits[exclude sky=region(source.reg)]'
              + '[bin time=::200]', opt='ltc1', outfile='lc_c7.fits',
              clobber='yes')
    deflare(infile='lc_c7.fits', outfile='outgood.gti', method='sigma',
            nsigma='2.5', plot='yes', save='lcsc_plot')
    dmcopy(infile=evt_file + '[@outgood.gti][ccd_id=7,sky=region(' + fov_file +
           '[ccd_id=7])]', outfile='evt2_df_ccd7.fits', clobber='yes')
    return


# Run the deflare script on the event file for a double source, excluding the
# source_1 region, which is the brightest
def deflareobs_double(evt_file, fov_file):
    dmcopy(infile=evt_file+'[energy=500:7000,ccd_id=7,' + 'sky=region('
           + fov_file + '[ccd_id=7])]', outfile='evt2_c7.fits', clobber='yes')
    print 'hit return to exit deflare plot'
    dmextract(infile='evt2_c7.fits[exclude sky=region(source_1.reg)]'
              + '[bin time=::200]', opt='ltc1', outfile='lc_c7.fits',
              clobber='yes')
    deflare(infile='lc_c7.fits', outfile='outgood.gti', method='sigma',
            nsigma='2.5', plot='yes', save='lcsc_plot')
    dmcopy(infile=evt_file + '[@outgood.gti][ccd_id=7,sky=region('
           + fov_file+'[ccd_id=7])]', outfile='evt2_df_ccd7.fits',
           clobber='yes')
    return
