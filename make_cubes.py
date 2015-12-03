import numpy as np 
import astropy.io.fits as pyfits
import astropy.wcs as pywcs
import numpy as np
import pylab
from scipy import ndimage
import glob
np.seterr(all='ignore')


def make_cube(filelist, outfilename, clobber=True, reference_image=None,
        reference_radec=None, reference_pix=None,rot90=False,cd1 = 0.95982,cd2 = 0.4472,cd1_2=0,cd2_1=0):

    """
    Create a data cube from a list of reduced TripleSpec spectra

    Parameters
    ----------
    filelist : list
        A list of the files to include.  Generate by, e.g.,
        filelist = glob.glob("*JHK.fits")
    outfilename : string
        The output FITS file name
    clobber : bool
        Overwrite the output file if it exists?
    reference_image : string
        A FITS file to use as the pointing standard for cross-correlation based
        pointing correction.  Makes use of agpy's cross_correlation_shift tool
    reference_radec : (float,float)
        The RA/Dec (CRVAL1/CRVAL2) of the image reference position.  If not
        specified, will be extracted from one of the input spectra (likely
        incorrectly)
    reference_pix : (float,float)
        The X,Y pixel (CRPIX1/CRPIX2) reference in the image.  If not
        specified, will default to the center.
    rot90 : bool
        Should the spectra be rotated by 90 degrees before stacking?
    """


    ### Find out the dimensions of the spectra
    ### fout is the HDU that the cube will be placed into
    fout = pyfits.open(filelist[0])
    spectral_axis_length = fout[0].header['NAXIS1']
    spatial_axis_length = fout[0].header['NAXIS2']
    number_of_files = len(filelist)

    ### Setup the cube and error cube
    cubeJHK = np.zeros([spatial_axis_length,spectral_axis_length,number_of_files])
    error_cubeJHK = np.zeros([spatial_axis_length,spectral_axis_length,number_of_files])
    
    ### Put each spectrum in the cube
    for ii,fn in enumerate(filelist): 
        d = pyfits.getdata(fn)
        cubeJHK[:,:,ii] = d[0,:,:]
        error_cubeJHK[:,:,ii] = d[1,:,:]
    print "cube shape: ",cubeJHK.shape

    ### Get reference RA/Dec    
    if reference_radec is not None:
        ra,dec = reference_radec
    else:
        fcen = pyfits.open(filelist[number_of_files/2])
        ra = fcen[0].header.get('RA')
        dec = fcen[0].header.get('DEC')
        ra,dec = coords.Position(ra+" "+dec).j2000()
    
    data = fout[0].data[0]
    error_data = fout[0].data[1] 
    data = cubeJHK = cubeJHK.swapaxes(0,1)
    error_data = error_cubeJHK = error_cubeJHK.swapaxes(0,1)

    ### If rotation is necessary, go.
    ### Also set spatial WCS parameters
    if not rot90:
        data = cubeJHK = data.swapaxes(1,2)[:,::-1,:]
        error_data = error_cubeJHK = error_data.swapaxes(1,2)[:,::-1,:]
        fout[0].header['CD1_1'] = -cd2/3600.
        fout[0].header['CD2_2'] =  cd1/3600.
        fout[0].header['CD1_2'] = cd1_2
        fout[0].header['CD2_1'] = cd2_1
    else:
        fout[0].header['CD1_1'] = -cd1/3600.
        fout[0].header['CD2_2'] =  cd2/3600.
        fout[0].header['CD1_2'] = cd2_1
        fout[0].header['CD2_1'] = cd1_2
    fout[0].header['NAXIS3'] = cubeJHK.shape[0]
    fout[0].header['NAXIS'] = 3
    print "fout.data shape: ",fout[0].data.shape

    if reference_pix is not None:
        xpix,ypix = reference_pix
    else:
        xpix = fout[0].data[0].shape[2]/2.+1
        ypix = fout[0].data[0].shape[1]


    #### Set more WCS parameters
    fout[0].header['CRVAL3'] = fout[0].header['CRVAL1']
    fout[0].header['CRPIX3'] = fout[0].header['CRPIX1']
    fout[0].header['CD3_3'] = fout[0].header['CDELT1']
    fout[0].header['CRVAL1'] = ra
    fout[0].header['CRVAL2'] = dec
    fout[0].header['CRPIX1'] = xpix
    fout[0].header['CRPIX2'] = ypix

    fout[0].header['CTYPE1'] = 'RA---TAN'
    fout[0].header['CTYPE2'] = 'DEC--TAN'
    fout[0].header['CTYPE3'] = 'LINEAR'
    del fout[0].header['CDELT1']
    del fout[0].header['CDELT2']
    

    #### Write cube to file
    outfilename_error = outfilename.replace(".fits","_error.fits")
    hdu_error = pyfits.PrimaryHDU(data=error_cubeJHK,header=fout[0].header)
    hdu_error.writeto(outfilename_error,clobber=clobber,output_verify='fix')
    fout[0].data = cubeJHK
    fout.writeto(outfilename,clobber=clobber,output_verify='fix')


if __name__ == "__main__":

#############
# Sample UT121124 slit scan
# The WCS solution parameters were found using IRAF ccmap, cctran, and ccsetwcs
    make_cube(sorted(glob.glob("6BNeast1*JHK_flat_backsub_airsub_cal_rvcorrect.fits")),
              "6BNeast1_UT121124_cube.fits",
              reference_radec=(83.8193901686863,-5.37339394961422),
              reference_pix=(53.7435709041544,40.2495441274706),
              cd2=.426,cd1=.505,cd1_2=-2.6616924140915e-6,cd2_1=3.52878723888732e-8,
              rot90=False)
##############


