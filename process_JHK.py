import numpy as np
import astropy.io.fits as pyfits
import pyspeckit
import glob
import astropy.stats.funcs
import pylab
import os
import sys

numcores=60

#########
## Read in the spectral-axis image coordinates for airglow lines and lines of interest 

coords_airglow_to_mask = np.loadtxt('airglow_lines_to_mask.txt',dtype='int') - 1 # shift IRAF image 
                                                                                 # coords to Python index coords
coords_airglow_to_subtract = np.loadtxt('airglow_lines_to_subtract.txt',dtype='int') - 1
coords_LOI = np.loadtxt('lines_of_interest_pixelcoords.txt',dtype='int') - 1


######### 
## Read in the telluric correction spectrum and define the correction to shift the spectrum into
                                                 ## the Kinematic Local Standard of Rest
calibrator_data = pyfits.getdata('HD37547_UT121124_xtellcor_tellspec_merged.fits')

VLSR = -9.03 ## km/s for Nov 24 2012 looking at Orion BN/KL


######### 
## Read-in the list of JHK spectra to process
filelist = sorted(glob.glob('Orion_BNKL_*_JHK.fits'))


######### 
## Begin for loop that cycles through each JHK spectrum in the list of JHK files
for h in range(len(filelist)):

    #####
    ## Open file and load data and header
    hdu = pyfits.open(filelist[h])
    data = hdu[0].data
    header = hdu[0].header

    #####
    ## copy the data for manipulation and setup the error array for error propagation along the way
    JHK_data = data.copy()
    error_data_total = data.copy()
    error_JHK_squared = data.copy()
    error_JHK_squared = error_JHK_squared/3.5 + (18./3.5/np.sqrt(header['NFOWLER']))**2
    
    ######### 
    ## read-in the normalized flat and flat-field the JHK spectrum
    BrQtz = pyfits.getdata('BrQtz_JHKnorm.fits')        
    JHK_data_flat = JHK_data / BrQtz[0]
    error_JHK_data_flat = np.sqrt( error_JHK_squared / BrQtz[0]**2 + BrQtz[1]**2 * (JHK_data / BrQtz[0]**2)**2 )
    
    #####
    ## Crop off the bottom and top spatial-direction pixels because they suffer from aliasing
    JHK_data_flat_cropped = JHK_data_flat[3:107,:]
    error_JHK_data_flat_cropped = error_JHK_data_flat[3:107,:]
    header['NAXIS2'] = 104
    
    #####
    ## Write the intermediate reduction product (flat-fielded and cropped) to file
    JHK_data_flat_cropped_witherror = np.array((JHK_data_flat_cropped,error_JHK_data_flat_cropped))
    hdu_flat = pyfits.PrimaryHDU(data=JHK_data_flat_cropped_witherror,header=header)
    filename = filelist[h].replace(".fits","_flat.fits")
    hdu_flat.writeto(filename,clobber=True)
 
    ############# 
    ## Begin masking
    JHK_data_flat_masked = JHK_data_flat_cropped.copy()
          
    #####
    ## Mask out rows with stars/continuum -- Don't judge me, I was young
    column_sum = np.sum(JHK_data_flat_masked[:,3891:3908],axis=1)
    stddev_column3898 = np.std(column_sum)
    if stddev_column3898 > 200.:
        masked_column3898 = astropy.stats.funcs.sigma_clip(column_sum,sig=1,iters=2)
        mask_stars = np.ma.getmask(masked_column3898)
        mask_index = np.where(mask_stars == True)
        mask_index_array = mask_index[0]

        numbers_to_append = ()
        for g in range(len(mask_index_array)):
            val = mask_index_array[g]
            if val == 0:
                numbers_to_append += (val,val+1,val+2)
            elif val == 1:
                numbers_to_append += (val-1,val,val+1,val+2)
            elif val == 102:
                numbers_to_append += (val-2,val-1,val,val+1)
            elif val == 103:
                numbers_to_append += (val-2,val-1,val)
            else:
                numbers_to_append += (val-2,val-1,val,val+1,val+2)

        numbers_to_append_array = np.array(numbers_to_append)
        
        JHK_data_flat_masked[numbers_to_append_array,:] = np.nan
        
    
    
        #### 
        ## Interpolate over NaNs from rows with stars
        for i in range(len(JHK_data_flat_masked[0,:])):
            data_column = JHK_data_flat_masked[:,i]
            mask = np.isnan(data_column)
            data_column[mask] = np.interp(np.flatnonzero(mask),np.flatnonzero(~mask),data_column[~mask])


    ######
    ## Remember for later where the zero-value region was in the JHK spectrum       
    index = np.where(JHK_data_flat_masked == 0)  
     
    ###### 
    ## Mask out airglow lines and lines of interest and interpolate over them  
    JHK_data_flat_masked[:,coords_airglow_to_mask_adjusted] = np.nan
    JHK_data_flat_masked[:,coords_LOI_adjusted] = np.nan
    for i in range(len(JHK_data_flat_masked[:,0])):
        data_row = JHK_data_flat_masked[i,:]
        mask = np.isnan(data_row)
        data_row[mask] = np.interp(np.flatnonzero(mask),np.flatnonzero(~mask),data_row[~mask])
   
    ##### 
    ## Set zero-value region back to 0
    JHK_data_flat_masked[index] = 0
        
    #####
    ## Write intermediate reduction product (masked and interpolated) to file
    hdu_flat_masked = pyfits.PrimaryHDU(data=JHK_data_flat_masked,header=header)
    filename = filelist[h].replace(".fits","_flat_masked.fits")
    hdu_flat_masked.writeto(filename,clobber=True)
    
    ############ 
    ## Fit the background -- this is ugly
    section_boundaries = [1,620,621,976,977,1230,1231,1750,1751,2512,2600,3629,3630,4647]            
    background_fit_2d = JHK_data_flat_masked.copy()
    error_background_fit_2d = JHK_data_flat_masked.copy()

    def fit_a_section(indices):
        i = indices[0]
        j = indices[1]
        range_lower,range_upper = section_boundaries[2*i],section_boundaries[2*i+1]
        
        if pylab.mod(j,10)==0:
            #check load
            m,n,o=os.getloadavg()
            if m>20:
                sys.exit()
        
        if i == 3:
            
            if j < 46:
                range_lower,range_upper = 1231,1750
            else:
                range_lower,range_upper = 1231,1749
                        
        if i == 4:
            
            if j < 46:
                range_lower,range_upper = 1751,2512
            else:
                range_lower,range_upper = 1750,2512
                        
        if i == 5:
            
            if j < 24:
                range_lower,range_upper = 2600,3630
            if j < 65 & j >= 24:
                range_lower,range_upper = 2600,3629
            if j > 65:
                range_lower,range_upper = 2600,3628
                    
        if i == 6:
            
            if j < 24:
                range_lower,range_upper = 3631,4647
            if j < 65 & j >= 24:
                range_lower,range_upper = 3630,4647
            if j > 65:
                range_lower,range_upper = 3629,4647
                
        section_background = JHK_data_flat_masked[j,range_lower-1:range_upper-1]
        wave_section_background = np.arange(len(JHK_data_flat_masked[0,range_lower-1:range_upper-1]))

        #####
        ## Setting up pyspeckit to fit a polynomial to each pixel row
        sp = pyspeckit.Spectrum(data=section_background,xarr=wave_section_background,header=pyfits.Header())
        sp.specfit.Registry.add_fitter('poly',pyspeckit.spectrum.models.polynomial_continuum.poly_fitter(order=2),
                                       3,multisingle='multi')
        #####
        ## Fitting
        sp.specfit(fittype='poly',guesses=[0,0,0], verbose=False)

        #####
        ## Results of fitting
        fit = sp.specfit.parinfo.values[0]*wave_section_background**2 + \
              sp.specfit.parinfo.values[1]*wave_section_background + sp.specfit.parinfo.values[2]
        error_fit_squared = (sp.specfit.parinfo.errors[0]*wave_section_background**2)**2 + \
                            (sp.specfit.parinfo.errors[1]*wave_section_background)**2 + \
                            sp.specfit.parinfo.errors[2]**2
        fit_witherrorsquared = np.array((fit,error_fit_squared))
        
        return fit_witherrorsquared


    #####
    ## Setup the indices for fitting        
    indices = []    ## Yuck
    for i in range(7):
        for j in range(104):
            indices.append([i,j])
    
    #####
    ## Run the fits in parallel            
    result = pyspeckit.parallel_map(fit_a_section,indices,numcores)

    #####
    ## Split up result to keep track of the fits
    sequence = np.array_split(indices,numcores)
    
    for k in range(len(sequence)):
      for y in range(len(sequence[k])):
        
        i = sequence[k][y][0]
        j = sequence[k][y][1]
        
        range_lower,range_upper = section_boundaries[2*i],section_boundaries[2*i+1]
        
        if i == 3:
            
            if j < 46:
                range_lower,range_upper = 1231,1750
            else:
                range_lower,range_upper = 1231,1749
                        
        if i == 4:
            
            if j < 46:
                range_lower,range_upper = 1751,2512
            else:
                range_lower,range_upper = 1750,2512
                        
        if i == 5:
            
            if j < 24:
                range_lower,range_upper = 2600,3630
            if j < 65 & j >= 24:
                range_lower,range_upper = 2600,3629
            if j > 65:
                range_lower,range_upper = 2600,3628
                    
        if i == 6:
            
            if j < 24:
                range_lower,range_upper = 3631,4647
            if j < 65 & j >= 24:
                range_lower,range_upper = 3630,4647
            if j > 65:
                range_lower,range_upper = 3629,4647
        
             
        background_fit_2d[j,range_lower-1:range_upper-1] = result[k][y][0]
        error_background_fit_2d[j,range_lower-1:range_upper-1] = np.sqrt(result[k][y][1])
    
    background_fit_2d_witherror = np.array((background_fit_2d,error_background_fit_2d))
    
    #####
    ## Save the background fits
    hdu_backfit = pyfits.PrimaryHDU(data=background_fit_2d_witherror,header=header)
    filename = filelist[h].replace(".fits","_flat_masked_backfit.fits")
    hdu_backfit.writeto(filename,clobber=True)    
    
    ########## 
    ## Subtract the background and write the background subtracted spectrum to file
    JHK_data_flat_backsub = JHK_data_flat_cropped - background_fit_2d
    error_JHK_data_flat_backsub = np.sqrt( error_JHK_data_flat_cropped**2 + error_background_fit_2d**2 )    
    JHK_data_flat_backsub_witherror = np.array((JHK_data_flat_backsub,error_JHK_data_flat_backsub))
    
    hdu_backsub = pyfits.PrimaryHDU(data=JHK_data_flat_backsub_witherror,header=header)
    filename = filelist[h].replace(".fits","_flat_backsub.fits")
    hdu_backsub.writeto(filename,clobber=True)
    
    #########
    ## Subtract airglow lines (not all airglow lines, just the ones 
    ## that don't overlap with lines of interest)

    ####
    ## Create indices for airglow lines (that are not blended with lines of interest) to be removed.            
    overlap_indices = np.in1d(coords_airglow_to_subtract_adjusted,coords_LOI_adjusted)    
    coords_airglow_adjusted_to_subtract = coords_airglow_to_subtract_adjusted[~overlap_indices]
    for i in range(len(coords_airglow_adjusted_to_subtract)):
        index = coords_airglow_adjusted_to_subtract[i]
        if stddev_column3898 > 200.:
            line = JHK_data_flat_backsub[:,index]
            median_to_subtract = np.median(line[~mask_stars])
        else:
            median_to_subtract = np.median(JHK_data_flat_backsub[:,index])
        JHK_data_flat_backsub[:,index] -= median_to_subtract

    ####
    ## Write airglow subtracted spectrum to file    
    JHK_data_flat_backsub_witherror = np.array((JHK_data_flat_backsub,error_JHK_data_flat_backsub))
    hdu_airsub = pyfits.PrimaryHDU(data=JHK_data_flat_backsub_witherror,header=header)
    filename = filelist[h].replace(".fits","_flat_backsub_airsub.fits")
    hdu_airsub.writeto(filename,clobber=True)
    
    
    ################    
    ## Flux calibration

    ###
    ## Create the wavelength array from the spectrum
    wavelength_array = (np.arange(header['NAXIS1']) - (header['CRPIX1'] - 1.) ) * header['CD1_1'] + header['CRVAL1']

    ###
    ## Interpolate the missing flux values so the wavelength array from the cube and that of the calibrator agree
    calibrator_interp_fluxvalues = np.interp(wavelength_array, calibrator_data[0,:]*1e4, calibrator_data[1,:])
    error_calibrator_interp_fluxvalues = np.interp(wavelength_array, calibrator_data[0,:]*1e4, calibrator_data[2,:])

    ###
    ## Multiply the data by the calibrator
    int_time = header['EXPTIME']
    JHK_data_flat_backsub_airsub_cal = JHK_data_flat_backsub.copy()
    error_JHK_data_flat_backsub_airsub_cal = JHK_data_flat_backsub.copy()    
    for i in range(len(JHK_data_flat_backsub[:,0])):
        JHK_data_flat_backsub_airsub_cal[i,:] = JHK_data_flat_backsub[i,:] * /
                                                calibrator_interp_fluxvalues / int_time
        error_JHK_data_flat_backsub_airsub_cal[i,:] = np.sqrt( (error_JHK_data_flat_backsub[i,:]* /
                                                      calibrator_interp_fluxvalues / int_time )**2 + /
                                                      (error_calibrator_interp_fluxvalues* /
                                                      JHK_data_flat_backsub[i,:]/int_time)**2 )
    

    ###
    ## Write flux-calibrated spectrum to file
    JHK_data_flat_backsub_airsub_cal_witherror = np.array((JHK_data_flat_backsub_airsub_cal,
                                                           error_JHK_data_flat_backsub_airsub_cal))    
    hdu_cal = pyfits.PrimaryHDU(data=JHK_data_flat_backsub_airsub_cal_witherror,header=header)
    filename = filelist[h].replace(".fits","_flat_backsub_airsub_cal.fits")
    hdu_cal.writeto(filename,clobber=True)    
    
    ######
    ## Shift the wavelength array into the Kinematic Local Standard of Rest    
    new_wavelength_array = wavelength_array * (1 + VLSR/3e5)
    JHK_data_flat_backsub_airsub_cal_rvcorrect = JHK_data_flat_backsub_airsub_cal.copy()
    error_JHK_data_flat_backsub_airsub_cal_rvcorrect = error_JHK_data_flat_backsub_airsub_cal.copy()
    
    for p in range(len(JHK_data_flat_backsub_airsub_cal[:,0])):
        interp_fluxvalues = np.interp(wavelength_array,new_wavelength_array,
                                      JHK_data_flat_backsub_airsub_cal[p,:])
        error_interp_fluxvalues = np.interp(wavelength_array,new_wavelength_array,
                                            error_JHK_data_flat_backsub_airsub_cal[p,:])
        JHK_data_flat_backsub_airsub_cal_rvcorrect[p,:] = interp_fluxvalues
        error_JHK_data_flat_backsub_airsub_cal_rvcorrect[p,:] = error_interp_fluxvalues
    


    ##########
    ## Write the FINAL spectrum (RV corrected) to file
    JHK_data_flat_backsub_airsub_cal_rvcorrect_witherror = /
                    np.array((JHK_data_flat_backsub_airsub_cal_rvcorrect,
                    error_JHK_data_flat_backsub_airsub_cal_rvcorrect))    
    hdu_rvcorrect = pyfits.PrimaryHDU(data=JHK_data_flat_backsub_airsub_cal_rvcorrect_witherror,header=header)
    filename = filelist[h].replace(".fits","_flat_backsub_airsub_cal_rvcorrect.fits")
    hdu_rvcorrect.writeto(filename,clobber=True)

