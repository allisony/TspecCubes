import astropy.io.fits as pyfits
import astropy.wcs as wcs
import numpy as np
import pyspeckit
from sdpy import makecube
import glob
from skimage.filter.rank import median as median_filter
from skimage.morphology import disk as sk_disk
import os
import astropy.units as u

#### Define prefix for the combined cube
pfx = 'Orion_TSPEC_jband'

#### Prefilter?
prefilter = False
if prefilter:
    cubefilename = pfx+"_cube_filtered.fits"
    nhitsfilename = pfx+"_cube_filtered_nhits.fits"
else:
    cubefilename = pfx+"_cube.fits"
    nhitsfilename = pfx+"_cube_nhits.fits"

####### For the 3 bands, these are the spectral axis WCS values to use
##K band: naxis3=2042, crval3=21776,cd3=2.8826265335083
##H band: naxis3=1511, crval3=16380,cd3=2.8826265335083
##J band: naxis3=998, crval3=12764, cd3=2.8826265335083

####### Generate the header for the combined cube
makecube.generate_header(83.809198, -5.3643852, naxis1=200, naxis2=250,
                         naxis3=998, coordsys='radec', ctype3='WAVE',
                         crval3=12764,
                         cd3=2.8826265335083,
                         cunit3='angstrom',
                         pixsize=1, bunit='erg/s/cm2/AA',
                         output_cubeheader='cube.hdr',
                         output_flatheader='flat.hdr',clobber=True)

if True: 
    cubeheader = pyfits.Header.fromtextfile('cube.hdr',endcard=False)
    cubewcs = wcs.WCS(cubeheader)

    cubeshape = [cubeheader['NAXIS%i' % i] for i in (3,2,1)]
    print "Making blank cube with size ",cubeshape
    blankcube = np.zeros(cubeshape,dtype='float32')

    cubefitsfile = pyfits.PrimaryHDU(data=blankcube, header=cubeheader)
    print "Writing blank cube to disk"
    cubefitsfile.writeto(cubefilename, clobber=True)

if True: 
    flatheader = pyfits.Header.fromtextfile('flat.hdr',endcard=False)
    blanknhits = np.empty([flatheader['NAXIS%i' % i] for i in (2,1)])
    nhitsfitsfile = pyfits.PrimaryHDU(data=blanknhits, header=flatheader)
    nhitsfitsfile.writeto(nhitsfilename, clobber=True)


### Filter function
def filter(fn, suffix='_filtered', clobber=True):
    fitsfile = pyfits.open(fn)
    data = fitsfile[0].data
    filtlist = [(data[ii,:,:]-np.median(data[ii,:,:],axis=0)[np.newaxis,:])
                for ii in xrange(data.shape[0])]
    filtlist = [im-im.min() for im in filtlist]
    filtered = np.array(filtlist)
    fitsfile[0].data = filtered
    outfn = fn.replace('.fits',suffix+".fits")
    if suffix not in outfn:
        raise ValueError("Not a FITS file?")
    else:
        fitsfile.writeto(outfn, clobber=clobber)
    return outfn

##### Combine the cubes
for fn in glob.glob("*cube.fits"):

    print "Working on file %s" % fn
    if pfx in fn:
        print "Skipping file ",fn
        continue

    hdr = pyfits.getheader(fn)
    W = wcs.WCS(pyspeckit.cubes.flatten_header(hdr))
    yinds,xinds = np.indices([hdr['NAXIS2'], hdr['NAXIS1']])
    wavel = hdr['CRVAL3'] + hdr['CD3_3']*(np.arange(hdr['NAXIS3'])+1-hdr['CRPIX3'])

    def data_iterator(data, **kwargs):
        for x,y in zip(xinds.ravel(),yinds.ravel()):
            yield data[:,y,x]

    def coord_iterator(data, **kwargs):
        for x,y in zip(xinds.ravel(),yinds.ravel()):
            ra,dec = W.wcs_pix2world(x,y,0)
            yield ra,dec

    def velo_iterator(data,linefreq=None, **kwargs):
        for x,y in zip(xinds.ravel(),yinds.ravel()):
            yield wavel

    if prefilter:
        fn = filter(fn)

    # this is where all the magic happens
    # add_with_kernel adds each pixel as a gaussian instead of just doing
    # nearest-neighbor interpolation
    # this will avoid empty stripes
    makecube.add_file_to_cube(fn,
                              cubefilename,
                              flatheader='flat.hdr',
                              cubeheader='cube.hdr', 
                              nhits=nhitsfilename,
                              data_iterator=data_iterator,
                              coord_iterator=coord_iterator,
                              velo_iterator=velo_iterator,
                              allow_smooth=False,
                              progressbar=True,
                              add_with_kernel=True,
                              kernel_fwhm=1.0/3600.,
                              default_unit=u.AA,
                              make_continuum=False)
                              #debug=1)

    print "Completed file %s" % fn
