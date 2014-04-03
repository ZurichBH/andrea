import matplotlib
matplotlib.use('Agg')
import aplpy
import os
import pyfits
import pywcs
import astropy.cosmology
import astropy.units as u
from numpy import *


# Extract the coordinate stored in the region file, double is a integer and
# can have values: {0,1,2}, where 0 is if the source is single, 1 for the
# first source in a double source and 2 for the second
def coord_source(ID, double):
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/'
             + 'Data/'+ID+'/repro/')
    if double == 0:
        with open('source.reg', 'r') as f:
            lines = f.read()
    elif double == 1:
        with open('source_1.reg', 'r') as f:
            lines = f.read()
    elif double == 2:
        with open('source_2.reg', 'r') as f:
            lines = f.read()
    ans = lines.find('circle')
    ans1 = ans+len('circle(')
    if lines[ans1+8] == ',':
        xpix = float(lines[ans1:ans1+8])
        if lines[ans1+17] == ',':
            ypix = float(lines[ans1+9:ans1+17])
        else:
            ypix = float(lines[ans1+9:ans1+18])
    else:
        xpix = float(lines[ans1:ans1+9])
        if lines[ans1+18] == ',':
            ypix = float(lines[ans1+10:ans1+18])
        else:
            ypix = float(lines[ans1+10:ans1+19])
    return xpix,  ypix


# Converts phyisical coordinates into wcs, using CHANDRA header for the
# reference point
def convert_pix2sky(ID, xpix, ypix):
    hdulist = pyfits.open('/Users/andyscanzio/Documents/ETH/Semestri/'
                          + 'FS2014/Project/Data/'+ID+'/repro/acisf'
                          + ID+'_repro_evt2.fits')
    wcs = pywcs.WCS(hdulist[1].header)
    tcrpx11 = hdulist[1].header['tcrpx11']
    tcrpx12 = hdulist[1].header['tcrpx12']
    tcdlt11 = hdulist[1].header['tcdlt11']
    tcdlt12 = hdulist[1].header['tcdlt12']
    tcrvl11 = hdulist[1].header['tcrvl11']
    tcrvl12 = hdulist[1].header['tcrvl12']
    tctyp11 = hdulist[1].header['tctyp11']
    tctyp12 = hdulist[1].header['tctyp12']
    wcs.wcs.crpix = [tcrpx11, tcrpx12]
    wcs.wcs.cdelt = array([tcdlt11, tcdlt12])
    wcs.wcs.crval = [tcrvl11, tcrvl12]
    wcs.wcs.ctype = [tctyp11, tctyp12]
    radec = wcs.wcs_pix2sky([[xpix, ypix]], 0)
    ra, dec = radec[0][0], radec[0][1]
    return ra, dec


# Plot the x-ray image, including the source, extended source and backgorund
# regions, the radius of the image is set so that the two galaxies are in the
# frame (with min = 65.0 and max = 155.0), double is used in the same way as
# in coord_source
def plot_xray(ID, double, z, dist, name):
    print ID
    gc = aplpy.FITSFigure('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                          + 'Project/Data/'+ID+'/repro/aplpy_'+ID
                          + '_0.5-8_flux.img')
    gc.show_colorscale(cmap='gist_heat', vmin=0, vmax=1.0e-07,
                       stretch='arcsinh', smooth=1)
    if double == 0:
        gc.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                        + 'Project/Data/'+ID+'/repro/source.reg')
        gc.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                        + 'Project/Data/'+ID+'/repro/source_ext.reg')
    elif double == 1:
        gc.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                        + 'Project/Data/'+ID+'/repro/source_1.reg')
        gc.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                        + 'Project/Data/'+ID+'/repro/source_1_ext.reg')
    elif double == 2:
        gc.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                        + 'Project/Data/'+ID+'/repro/source_2.reg')
        gc.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                        + 'Project/Data/'+ID+'/repro/source_2_ext.reg')
    gc.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/'
                    + 'Data/'+ID+'/repro/back.reg')

    radius = dist*u.kpc*astropy.cosmology.arcsec_per_kpc_comoving(z)/u.arcsec
    radius = radius.value
    print radius
    if radius < 60.0:
        radius = 60.0
    elif radius > 150.0:
        radius = 150.0
    radius += 5.0

    scale = 10.0*u.kpc/astropy.cosmology.arcsec_per_kpc_comoving(z)
    scale = scale.value
    gc.add_label(0.20, 0.1, 'ObsID: '+ID+'\n'+name, relative=True,
                 weight='heavy', size='x-large', color='gray')

    xpix, ypix = coord_source(ID, double)
    ra, dec = convert_pix2sky(ID, xpix, ypix)

    gc.recenter(ra, dec, radius=radius/3600.0)
    gc.tick_labels.set_xformat('ddmmss')
    gc.tick_labels.set_yformat('ddmmss')
    gc.add_scalebar(10.0/3600.0, linewidth=3.0)
    gc.add_colorbar()
    gc.scalebar.set_label('10.0"/'+str(scale)[0:4]+' kpc')
    gc.scalebar.set_font(size='x-large', weight='heavy')
    gc.scalebar.set_color('gray')
    if double == 0:
        gc.save('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                + 'Project/Images/'+ID+'_xray.eps')
    elif double == 1:
        gc.save('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                + 'Project/Images/'+ID+'_1_xray.eps')
    elif double == 2:
        gc.save('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                + 'Project/Images/'+ID+'_2_xray.eps')
    gc.close()


# Plot the optical image with source and extended source regions, using the
# same settings as plot_xray
def plot_optical(ID, double, z, dist, name):
    print ID
    gc = aplpy.FITSFigure('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                          + 'Project/Optical/'+ID+'_optical.fits')
    gc.show_grayscale()

    gc.add_label(0.20, 0.1, 'ObsID: '+ID+'\n'+name, relative=True,
                 weight='heavy', size='x-large', color='gray')

    radius = dist*u.kpc*astropy.cosmology.arcsec_per_kpc_comoving(z)/u.arcsec
    radius = radius.value
    print radius
    if radius < 60.0:
        radius = 60.0
    elif radius > 150.0:
        radius = 150.0
    radius += 5.0

    scale = 10.0*u.kpc/astropy.cosmology.arcsec_per_kpc_comoving(z)
    scale = scale.value

    xpix, ypix = coord_source(ID, double)
    ra, dec = convert_pix2sky(ID, xpix, ypix)

    gc.show_circles(ra, dec, 1.5/3600.0, color='green', linewidth=2.0)
    ri = 2.0*u.kpc*astropy.cosmology.arcsec_per_kpc_comoving(z)/u.arcsec
    ri = ri.value
    ro = 4.0*u.kpc*astropy.cosmology.arcsec_per_kpc_comoving(z)/u.arcsec
    ro = ro.value
    gc.show_circles(ra, dec, ri/3600.0, color='cyan', linewidth=2.0)
    gc.show_circles(ra, dec, ro/3600.0, color='cyan', linewidth=2.0)

    gc.recenter(ra, dec, radius=radius/3600.0)
    gc.tick_labels.set_xformat('ddmmss')
    gc.tick_labels.set_yformat('ddmmss')
    gc.add_scalebar(10.0/3600.0, linewidth=3.0)

    gc.scalebar.set_label('10.0"/'+str(scale)[0:4]+' kpc')
    gc.scalebar.set_font(size='x-large', weight='heavy')
    gc.scalebar.set_color('gray')
    if double == 0:
        gc.save('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/'
                + 'Images/'+ID+'_optical.eps')
    elif double == 1:
        gc.save('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/'
                + 'Images/'+ID+'_1_optical.eps')
    elif double == 2:
        gc.save('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/'
                + 'Images/'+ID+'_2_optical.eps')
    gc.close()


# Calculate the physical area in kpc of the extended source region
def area_spec(ID, z):
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Data/'
             + ID[0:5]+'/repro/')
    area = pyfits.getval(ID+'_spec_ext.pi', 'backscal', ext=1)*8192**2
    area = area*0.25*u.arcsec**2
    area = area/(astropy.cosmology.arcsec_per_kpc_comoving(z)**2)
    area = float(area/u.kpc**2)
    return area
