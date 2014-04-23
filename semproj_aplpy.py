import matplotlib
matplotlib.use('Agg')
import aplpy
import os
import pyfits
import pywcs
import astropy.cosmology
import astropy.units as u
from numpy import *
from matplotlib.pyplot import *


def opt_regions(gc, ID, double, z):
    xpix, ypix = coord_source(ID, double)
    ra, dec = convert_pix2sky(ID, xpix, ypix)

    ri = 2.0*u.kpc*astropy.cosmology.arcsec_per_kpc_comoving(z)/u.arcsec
    ri = ri.value
    ro = 4.0*u.kpc*astropy.cosmology.arcsec_per_kpc_comoving(z)/u.arcsec
    ro = ro.value
    if double == 2:
        gc.show_circles(ra, dec, 1.5/3600.0, color='blue', linewidth=2.0)
        gc.show_circles(ra, dec, ri/3600.0, color='magenta', linewidth=2.0)
        gc.show_circles(ra, dec, ro/3600.0, color='magenta', linewidth=2.0)
    else:
        gc.show_circles(ra, dec, 1.5/3600.0, color='green', linewidth=2.0)
        gc.show_circles(ra, dec, ri/3600.0, color='cyan', linewidth=2.0)
        gc.show_circles(ra, dec, ro/3600.0, color='cyan', linewidth=2.0)


def settings(gc, ID, name, z, ra, dec, dist):
    radius = dist*u.kpc*astropy.cosmology.arcsec_per_kpc_comoving(z)/u.arcsec
    radius = radius.value
    if radius < 15.0:
        radius = 15.0
    elif radius > 225.0:
        radius = 225.0
    radius += 15.0

    scale = 10.0*u.kpc/astropy.cosmology.arcsec_per_kpc_comoving(z)
    scale = scale.value

    gc.recenter(ra, dec, radius=radius/3600.0)
    gc.tick_labels.set_xformat('ddmmss')
    gc.tick_labels.set_yformat('ddmmss')
    gc.add_scalebar(10.0/3600.0, linewidth=3.0)
    gc.scalebar.set_label('10.0"/'+str(scale)[0:4]+' kpc')
    gc.scalebar.set_font(size='x-large', weight='heavy')
    gc.scalebar.set_color('gray')


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


def plot_images(ID, double, z, dist, name):
    print ID

    fig = figure(figsize=(20, 10))
    f1 = aplpy.FITSFigure('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                          + 'Project/Data/'+ID+'/repro/aplpy_'+ID
                          + '_0.5-8_flux.img', figure=fig, subplot=(1, 2, 1))
    f1.show_colorscale(cmap='gist_heat', vmin=0, pmax=99.9,
                       stretch='arcsinh', smooth=1)
    if double == 0:
        f1.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                        + 'Project/Data/'+ID+'/repro/source.reg')
        f1.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                        + 'Project/Data/'+ID+'/repro/source_ext.reg')
    elif double == 1:
        f1.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                        + 'Project/Data/'+ID+'/repro/source_1.reg')
        f1.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                        + 'Project/Data/'+ID+'/repro/source_1_ext.reg')
        f1.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                        + 'Project/Data/'+ID+'/repro/source_2.reg')
        f1.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                        + 'Project/Data/'+ID+'/repro/source_2_ext.reg')

    xpix, ypix = coord_source(ID, double)
    ra, dec = convert_pix2sky(ID, xpix, ypix)

    settings(f1, ID, name, z, ra, dec, dist)

    f2 = aplpy.FITSFigure('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                          + 'Project/Optical/'+ID+'_optical.fits', figure=fig,
                          subplot=(1, 2, 2))
    f2.show_grayscale(stretch='arcsinh', pmax=99.9)

    opt_regions(f2, ID, double, z)
    if double == 1:
        opt_regions(f2, ID, 2, z)

    settings(f2, ID, name, z, ra, dec, dist)

    suptitle(ID+' - '+name, size='xx-large', weight='heavy')
    fig.text(0.285, 0.87, 'X-ray', size='large', weight='bold')
    fig.text(0.715, 0.87, 'Optical', size='large', weight='bold')

    fig.savefig('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
                + 'Project/Images/'+ID+'_images.eps')


# Calculate the physical area in kpc of the extended source region
def area_spec(ID, z):
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Data/'
             + ID[0:5]+'/repro/')
    area = pyfits.getval(ID+'_spec_ext.pi', 'backscal', ext=1)*8192**2
    area = area*0.25*u.arcsec**2
    area = area/(astropy.cosmology.arcsec_per_kpc_comoving(z)**2)
    area = float(area/u.kpc**2)
    return area


def exptime(ID):
    hdulist = pyfits.open('/Users/andyscanzio/Documents/ETH/Semestri/'
                          + 'FS2014/Project/Data/'+ID+'/repro/acisf'
                          + ID+'_repro_evt2.fits')
    return hdulist[1].header['EXPTIME']
