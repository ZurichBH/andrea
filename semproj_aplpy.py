import matplotlib
matplotlib.use('Agg')
import aplpy
import os
import pyfits
import pywcs
import astropy.cosmology
import astropy.units as u


# Extract the physical coordinates from a region file, using the keyword "circle"
# as a starting point for where the coordinates begins.
def coord_source(ID):
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Data/'
             + ID+'/repro/')
    with open('source.reg', 'r') as f:
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


# Same as coord_source, but for the first source in a double source
def coord_source_1(ID):
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Data/'
             + ID + '/repro/')
    with open('source_1.reg', 'r') as f:
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


# Same as coord_source, but for the second source in a double source
def coord_source_2(ID):
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Data/'
             + ID + '/repro/')
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


def plot_xray(ID, z, dist, name):
    print ID
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Data/'
             + ID + '/repro/')
    gc = aplpy.FITSFigure('aplpy_'+ID+'_0.5-8_flux.img')
    gc.show_colorscale(cmap='gist_heat', vmin=0, vmax=1.0e-07,
                       stretch='arcsinh', smooth=1)
    gc.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/'
                    + 'Data/'+ID+'/repro/source.reg')
    gc.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/'
                    + 'Data/'+ID+'/repro/source_ext.reg')
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

    xpix, ypix = coord_source(ID)
    hdulist = pyfits.open('/Users/andyscanzio/Documents/ETH/Semestri/'
                          + 'FS2014/Project/Data/'+ID+'/repro/acisf'
                          + ID+'_repro_evt2.fits')
    wcs = pywcs.WCS(hdulist[1].header)
    print wcs.wcs.name
    wcs.wcs.print_contents()
    radec = wcs.wcs_pix2sky([[xpix, ypix]], 0)
    print radec
    ra, dec = radec[0][0], radec[0][1]
    print ra, dec

    gc.recenter(ra, dec, radius=radius/3600.0)
    gc.tick_labels.set_xformat('ddmmss')
    gc.tick_labels.set_yformat('ddmmss')
    gc.add_scalebar(10.0/3600.0)
    gc.add_colorbar()
    gc.scalebar.set_label('10.0"/'+str(scale)[0:4]+' kpc')
    gc.scalebar.set_font(size='x-large', weight='heavy')
    gc.scalebar.set_color('gray')
    gc.save('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/'
            + 'Project/Images/'+ID+'_xray.eps')
    gc.close()


def plot_xray_1(ID, z, dist, name):
    print ID
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Data/'
             + ID+'/repro/')
    gc = aplpy.FITSFigure('aplpy_'+ID+'_0.5-8_flux.img')
    gc.show_colorscale(cmap='gist_heat', vmin=0, vmax=1.0e-07,
                       stretch='arcsinh', smooth=1)
    gc.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/'
                    + 'Data/'+ID+'/repro/source_1.reg')
    gc.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/'
                    + 'Data/'+ID+'/repro/source_1_ext.reg')
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

    xpix, ypix = coord_source(ID)
    hdulist = pyfits.open('/Users/andyscanzio/Documents/ETH/Semestri/'
                          + 'FS2014/Project/Data/'+ID+'/repro/acisf'
                          + ID+'_repro_evt2.fits')
    wcs = pywcs.WCS(hdulist[0].header)
    radec = wcs.wcs_pix2sky([[xpix, ypix]], 0)
    ra, dec = radec[0][0], radec[0][1]

    gc.recenter(ra, dec, radius=radius/3600.0)
    gc.tick_labels.set_xformat('ddmmss')
    gc.tick_labels.set_yformat('ddmmss')
    gc.add_scalebar(10.0/3600.0)
    gc.add_colorbar()
    gc.scalebar.set_label('10.0"/'+str(scale)[0:4]+' kpc')
    gc.scalebar.set_font(size='x-large', weight='heavy')
    gc.scalebar.set_color('gray')
    gc.save('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Images/'
            + ID+'_1_xray.eps')
    gc.close()


def plot_xray_2(ID, z, dist, name):
    print ID
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Data/'
             + ID+'/repro/')
    gc = aplpy.FITSFigure('aplpy_'+ID+'_0.5-8_flux.img')
    gc.show_colorscale(cmap='gist_heat', vmin=0, vmax=1.0e-07,
                       stretch='arcsinh', smooth=1)
    gc.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/'
                    + 'Data/'+ID+'/repro/source_2.reg')
    gc.show_regions('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/'
                    + 'Data/'+ID+'/repro/source_2_ext.reg')
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

    xpix, ypix = coord_source(ID)
    hdulist = pyfits.open('/Users/andyscanzio/Documents/ETH/Semestri/'
                          + 'FS2014/Project/Data/'+ID+'/repro/acisf'
                          + ID+'_repro_evt2.fits')
    wcs = pywcs.WCS(hdulist[0].header)
    radec = wcs.wcs_pix2sky([[xpix, ypix]], 0)
    ra, dec = radec[0][0], radec[0][1]

    gc.recenter(ra, dec, radius=radius/3600.0)
    gc.tick_labels.set_xformat('ddmmss')
    gc.tick_labels.set_yformat('ddmmss')
    gc.add_scalebar(10.0/3600.0)
    gc.add_colorbar()
    gc.scalebar.set_label('10.0"/'+str(scale)[0:4]+' kpc')
    gc.scalebar.set_font(size='x-large', weight='heavy')
    gc.scalebar.set_color('gray')
    gc.save('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Images/'
            + ID+'_2_xray.eps')
    gc.close()


def plot_optical(ID, z, dist, name):
    print ID
    gc = aplpy.FITSFigure(ID+'_optical.fits')
    gc.show_grayscale()

    xpix, ypix = coord_source(ID)
    hdulist = pyfits.open('/Users/andyscanzio/Documents/ETH/Semestri/'
                          + 'FS2014/Project/Data/'+ID+'/repro/acisf'
                          + ID+'_repro_evt2.fits')
    wcs = pywcs.WCS(hdulist[0].header)
    radec = wcs.wcs_pix2sky([[xpix, ypix]], 0)
    ra, dec = radec[0][0], radec[0][1]

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

    gc.recenter(ra, dec, radius=radius/3600.0)
    gc.tick_labels.set_xformat('ddmmss')
    gc.tick_labels.set_yformat('ddmmss')
    gc.add_scalebar(10.0/3600.0)

    gc.scalebar.set_label('10.0"/'+str(scale)[0:4]+' kpc')
    gc.scalebar.set_font(size='x-large', weight='heavy')
    gc.scalebar.set_color('gray')
    gc.save('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Images/'
            + ID+'_optical.eps')
    gc.close()


def plot_optical_1(ID, z, dist, name):
    print ID
    gc = aplpy.FITSFigure(ID+'_optical.fits')
    gc.show_grayscale()

    xpix, ypix = coord_source(ID)
    hdulist = pyfits.open('/Users/andyscanzio/Documents/ETH/Semestri/'
                          + 'FS2014/Project/Data/'+ID+'/repro/acisf'
                          + ID+'_repro_evt2.fits')
    wcs = pywcs.WCS(hdulist[0].header)
    radec = wcs.wcs_pix2sky([[xpix, ypix]], 0)
    ra, dec = radec[0][0], radec[0][1]

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

    gc.recenter(ra, dec, radius=radius/3600.0)
    gc.tick_labels.set_xformat('ddmmss')
    gc.tick_labels.set_yformat('ddmmss')
    gc.add_scalebar(10.0/3600.0)

    gc.scalebar.set_label('10.0"/'+str(scale)[0:4]+' kpc')
    gc.scalebar.set_font(size='x-large', weight='heavy')
    gc.scalebar.set_color('gray')
    gc.save('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Images/'
            + ID+'_1_optical.eps')
    gc.close()


def plot_optical_2(ID, z, dist, name):
    print ID
    gc = aplpy.FITSFigure(ID+'_optical.fits')
    gc.show_grayscale()

    xpix, ypix = coord_source(ID)
    hdulist = pyfits.open('/Users/andyscanzio/Documents/ETH/Semestri/'
                          + 'FS2014/Project/Data/'+ID+'/repro/acisf'
                          + ID+'_repro_evt2.fits')
    wcs = pywcs.WCS(hdulist[0].header)
    radec = wcs.wcs_pix2sky([[xpix, ypix]], 0)
    ra, dec = radec[0][0], radec[0][1]

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

    gc.recenter(ra, dec, radius=radius/3600.0)
    gc.tick_labels.set_xformat('ddmmss')
    gc.tick_labels.set_yformat('ddmmss')
    gc.add_scalebar(10.0/3600.0)

    gc.scalebar.set_label('10.0"/'+str(scale)[0:4]+' kpc')
    gc.scalebar.set_font(size='x-large', weight='heavy')
    gc.scalebar.set_color('gray')
    gc.save('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Images/'
            + ID+'_2_optical.eps')
    gc.close()


def area_spec(ID, z):
    os.chdir('/Users/andyscanzio/Documents/ETH/Semestri/FS2014/Project/Data/'
             + ID[0:5]+'/repro/')
    area = pyfits.getval(ID+'_spec_ext.pi', 'backscal', ext=1)*8192**2
    area = area*0.25*u.arcsec**2
    area = area/(astropy.cosmology.arcsec_per_kpc_comoving(z)**2)
    area = float(area/u.kpc**2)
    return area
