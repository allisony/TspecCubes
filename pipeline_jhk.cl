procedure pipeline_jhk (filename,outname)
string filename
string outname
string procdir='./'             {prompt="Output directory for intermediate files"}
string rawdir='./'              {prompt="Directory in which the raw files are stored"}
string dark=''                  {prompt="Filename of dark (assumes current dir)"}
string skylines='skylines'      {prompt="Filename of fcskylines file to use (assumes current dir)"}
bool darksub=yes                {prompt="Do dark subtraction? (must specify dark!)"}
bool recrop=yes                 {prompt="Redo cropping / transforming?"}
bool justcrop=no                {prompt="Just crop (no transform?)"}
bool combine=yes                {prompt="Combine images?"}
bool background=yes             {prompt="Run background subtraction?"}
bool calibrate=no               {prompt="Calibrate?"}
bool redowavelength=no          {prompt="Run reidentify for wavelength solution?"}
string search="5"               {prompt="Search radius for reidentify (INDEF or 5 have been tried)"}
bool interactive_background=no  {prompt="Do background subtraction interactively?"}
bool calibrator=no              {prompt="Is the source a calibrator?"}
bool clobber=yes                {prompt="Overwrite existing files?"}
bool caleach=no                 {prompt="Calibrate EACH file? (requires pre-existing calibrator data)"}
bool noJ=no                     {prompt="Ignore J band?"}
string sample="*"               {prompt="Selection region for background task"}
real magJ=0.0                   {prompt="J Magnitude of calibrator star"}
real magH=0.0                   {prompt="H Magnitude of calibrator star"}
real magK=0.0                   {prompt="K Magnitude of calibrator star"}
real teff=10000                 {prompt="Effective temperature of calibrator (required for calibrator=yes)"}
struct *flist
struct *flist2
begin
    struct line
    string filenameJ
    string filenameH
    string filenameK
    string ds
    real KdeltaL,HdeltaL,JdeltaL,KdeltaP,HdeltaP,JdeltaP
    flist=filename
    flist2=filename
    delete ( filename+"J" )
    delete ( filename+"H" )
    delete ( filename+"K" )
    print ( "", > filename+"J" )
    print ( "", > filename+"H" )
    print ( "", > filename+"K" )

    onedspec
    twodspec
    apextract
    longslit

    if (recrop) {
    print ("*****Beginning recrop loop*****")
        while(fscan(flist,line)!=EOF) {
            fixpix ( rawdir+line , "bpm.fits" )
            hedit ( rawdir+line,"DISPAXIS",1,add+,ver-)
            if (darksub) { 
                imarith ( rawdir+line, '-', dark, procdir+line+'_darksub.fits' ) 
                ds = '_darksub' 
                print ( "Subtracted dark "+dark+" from "+line )
                # separate out the individual orders
                imcopy (procdir+line+"_darksub[*,710:833]", procdir+line+"_K.fits", verbose="yes")
                imcopy (procdir+line+"_darksub[*,532:665]", procdir+line+"_H.fits", verbose="yes")
                imcopy (procdir+line+"_darksub[*,350:537]", procdir+line+"_J.fits", verbose="yes")
            } else {
                imcopy (rawdir+line+"[*,710:833]", procdir+line+"_K.fits", verbose="yes")
                imcopy (rawdir+line+"[*,532:665]", procdir+line+"_H.fits", verbose="yes")
                imcopy (rawdir+line+"[*,350:537]", procdir+line+"_J.fits", verbose="yes")
            }

            filenameJ = line+"_J"
            filenameH = line+"_H"
            filenameK = line+"_K"

            # apply pre-calculated coordinate transforms
            if (justcrop==no) {
                # reidentify will fail if you specify a procdir.  I don't know if there's a workaround for this
                # I suspect not =(
                if (redowavelength==yes) {
                    reidentify ( "sample_id_J" , procdir+filenameJ,  coordlist="OHll.dat", nlost=10, verbose="yes", shift="INDEF", search=search, step=5, nsum=5, override="yes" )
                    reidentify ( "sample_id_H" , procdir+filenameH,  coordlist="OHll.dat", nlost=10, verbose="yes", shift="INDEF", search=search, step=5, nsum=5, override="yes" )
                    reidentify ( "sample_id_K" , procdir+filenameK,  coordlist="OHll.dat", nlost=10, verbose="yes", shift="INDEF", search=search, step=5, nsum=5, override="yes" )
                    fitcoords ( procdir+filenameJ, fitname=filenameJ+"_skylinesJ", combine=yes, xo=4, yo=3, interactive="no" )
                    fitcoords ( procdir+filenameH, fitname=filenameH+"_skylinesH", combine=yes, xo=4, yo=3, interactive="no" )
                    fitcoords ( procdir+filenameK, fitname=filenameK+"_skylinesK", combine=yes, xo=4, yo=3, interactive="no" )
                    transform ( input=procdir+filenameJ , output=procdir+filenameJ+"t", fitnames=filenameJ+"_skylinesJ,starsJ" )
                    transform ( input=procdir+filenameH , output=procdir+filenameH+"t", fitnames=filenameH+"_skylinesH,starsH" )
                    transform ( input=procdir+filenameK , output=procdir+filenameK+"t", fitnames=filenameK+"_skylinesK,starsK" )
                }
                else {
                    #transform ( input=procdir+filenameJ , output=procdir+filenameJ+"t", fitnames=skylines+"J,starsJ" )
                    #transform ( input=procdir+filenameH , output=procdir+filenameH+"t", fitnames=skylines+"H,starsH" )
                    #transform ( input=procdir+filenameK , output=procdir+filenameK+"t", fitnames=skylines+"K,starsK" )

                    print (procdir+filenameJ, procdir+filenameJ+"t", skylines+"J,starsJ")
                    transform ( input=procdir+filenameJ , output=procdir+filenameJ+"t", fitnames=skylines+"J,starsJ", x1=11324.6533203125, dx=2.88262653351, nx=1230, dy=1.0, database='database' )
                    transform ( input=procdir+filenameH , output=procdir+filenameH+"t", fitnames=skylines+"H,starsH", x1=14138.3642578125, dx=2.88262653351, nx=1536, dy=1.0, database='database')
                    transform ( input=procdir+filenameK , output=procdir+filenameK+"t", fitnames=skylines+"K,starsK", x1=18817.3847656, dx=2.88262653351, nx=2048, dy=1.0, database='database')
                }
                imcopy (procdir+filenameJ+"t[*,69:178]", procdir+filenameJ+"tc.fits", verbose="yes" )
                imcopy (procdir+filenameH+"t[*,21:130]", procdir+filenameH+"tc.fits", verbose="yes" )
                imcopy (procdir+filenameK+"t[*,16:125]",   procdir+filenameK+"tc.fits", verbose="yes" )
                

                hedit (procdir+filenameJ+"tc.fits", "CRVAL2", 1., verify=no, update=yes)
                hedit (procdir+filenameJ+"tc.fits", "CRPIX2", 1., verify=no, update=yes)
                hedit (procdir+filenameH+"tc.fits", "CRVAL2", 1., verify=no, update=yes)
                hedit (procdir+filenameH+"tc.fits", "CRPIX2", 1., verify=no, update=yes)
                hedit (procdir+filenameK+"tc.fits", "CRVAL2", 1., verify=no, update=yes)
                hedit (procdir+filenameK+"tc.fits", "CRPIX2", 1., verify=no, update=yes)

                print ( procdir+filenameJ+"tc", >> filename+"J" )
                print ( procdir+filenameH+"tc", >> filename+"H" )
                print ( procdir+filenameK+"tc", >> filename+"K" )
            }

        }
      print ("*****Completed recrop loop*****")
    }

    # calibrate using some sort of standard/sensfunc (I used blackbodies)
    if (caleach) {
    print ("*****Beginning caleach loop*****")
         while(fscan(flist2,line)!=EOF) {
            filenameJ = line+"_J"
            filenameH = line+"_H"
            filenameK = line+"_K"


            if (clobber) {
                if (noJ) {
                    delete ( procdir+filename+"_HK.fits"  )
                }
                else {
                    delete ( procdir+filename+"_JHK.fits"  )
                }
            }

            if (noJ) {
                imcombine ( procdir+filenameK+"tc,"+procdir+filenameH+"tc" , procdir+line+"_HK" , combine="sum" , offset="wcs" )
            }
            else {
                imcombine ( procdir+filenameK+"tc,"+procdir+filenameH+"tc,"+procdir+filenameJ+"tc" , procdir+line+"_JHK" , combine="sum" , offset="wcs" )
            }
            if (background) {
                if (clobber) {
                    if (noJ) {
                        delete ( procdir+line+"_HK_backsub.fits" )
                    }
                    else {
                        delete ( procdir+line+"_JHK_backsub.fits" )
                    }
                }
                if (noJ) {
                    background ( procdir+line+"_HK" , procdir+line+"_HK_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background, sample=sample )
                } else {
                    background ( procdir+line+"_JHK" , procdir+line+"_JHK_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background, sample=sample )
                }
                print ("Background subtraction complete!")
            }
        }
     print ("*****Completed caleach loop*****")
    }

    # combine J/H/K spectra
    if (combine) {
        if (! noJ) {
            imcombine ( "@"+filename+"J" , procdir+outname+"_J_combine" , combine="median" , scale="mode" )
        }
        imcombine ( "@"+filename+"H" , procdir+outname+"_H_combine" , combine="median" , scale="mode" )
        imcombine ( "@"+filename+"K" , procdir+outname+"_K_combine" , combine="median" , scale="mode" )
    }


    # calibrate using some sort of standard/sensfunc (I used blackbodies)
    if (! caleach) {
        if (calibrate) {
            calibrate ( procdir+outname+"_J_combine" , procdir+outname+"_J_cal" , sens="sensJ" , ignoreaps+)
            calibrate ( procdir+outname+"_H_combine" , procdir+outname+"_H_cal" , sens="sensH" , ignoreaps+)
            calibrate ( procdir+outname+"_K_combine" , procdir+outname+"_K_cal" , sens="sensK" , ignoreaps+)
            imgets ( procdir+outname+"_K_cal" , "CD1_1")
            KdeltaL = imgets.value
            imgets ( procdir+outname+"_J_cal" , "CD1_1")
            JdeltaL = imgets.value
            imgets ( procdir+outname+"_H_cal" , "CD1_1")
            HdeltaL = imgets.value
            imgets ( procdir+outname+"_K_cal" , "CD2_2")
            KdeltaP = imgets.value
            imgets ( procdir+outname+"_J_cal" , "CD2_2")
            JdeltaP = imgets.value
            imgets ( procdir+outname+"_H_cal" , "CD2_2")
            HdeltaP = imgets.value
            magnify ( procdir+outname+"_J_cal" , procdir+outname+"_J_cal_m" , JdeltaL/KdeltaL , JdeltaP/KdeltaP )
            magnify ( procdir+outname+"_H_cal" , procdir+outname+"_H_cal_m" , HdeltaL/KdeltaL , HdeltaP/KdeltaP )
            if (noJ) {
                imcombine ( procdir+outname+"_K_cal,"+procdir+outname+"_H_cal_m", procdir+outname+"_HK_cal" , combine="sum" , offset="wcs" )
                if (background) {
                    background ( procdir+outname+"_HK_cal" , procdir+outname+"_HK_cal_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background, sample=sample )
                }
            }
            else {
                imcombine ( procdir+outname+"_K_cal,"+procdir+outname+"_H_cal_m,"+procdir+outname+"_J_cal_m" , procdir+outname+"_JHK_cal" , combine="sum" , offset="wcs" )
                if (background) {
                    background ( procdir+outname+"_JHK_cal" , procdir+outname+"_JHK_cal_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background, sample=sample )
                }
            }
        }
    # background subtract to remove night sky lines
        else {
            if (background) {
                imgets ( procdir+outname+"_K_combine" , "CD1_1")
                KdeltaL = imgets.value
                imgets ( procdir+outname+"_J_combine" , "CD1_1")
                JdeltaL = imgets.value
                imgets ( procdir+outname+"_H_combine" , "CD1_1")
                HdeltaL = imgets.value
                imgets ( procdir+outname+"_K_combine" , "CD2_2")
                KdeltaP = imgets.value
                imgets ( procdir+outname+"_J_combine" , "CD2_2")
                JdeltaP = imgets.value
                imgets ( procdir+outname+"_H_combine" , "CD2_2")
                HdeltaP = imgets.value
                magnify ( procdir+outname+"_J_combine" , procdir+outname+"_J_combine_m" , JdeltaL/KdeltaL , JdeltaP/KdeltaP )
                magnify ( procdir+outname+"_H_combine" , procdir+outname+"_H_combine_m" , HdeltaL/KdeltaL , HdeltaP/KdeltaP )
                if (calibrator) {
                    if (background) {
                        background ( procdir+outname+"_J_combine" , procdir+outname+"_J_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background, sample=sample )
                        background ( procdir+outname+"_H_combine" , procdir+outname+"_H_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background, sample=sample )
                        background ( procdir+outname+"_K_combine" , procdir+outname+"_K_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background, sample=sample )
                    }
                    apall ( procdir+outname+"_J_backsub" , interactive="no" )
                    apall ( procdir+outname+"_H_backsub" , interactive="no" )
                    apall ( procdir+outname+"_K_backsub" , interactive="no" )
                    standard ( procdir+outname+"_J_backsub.ms" , output="stdJ" , star_name="J" , caldir="onedstds$blackbody/" , mag=magJ, teff=teff , magband="J" )
                    standard ( procdir+outname+"_H_backsub.ms" , output="stdH" , star_name="H" , caldir="onedstds$blackbody/" , mag=magH, teff=teff , magband="H" )
                    standard ( procdir+outname+"_K_backsub.ms" , output="stdK" , star_name="K" , caldir="onedstds$blackbody/" , mag=magK, teff=teff , magband="K" )
                    sensfunc ( "stdJ" , "sensJ" , answer="YES" )
                    sensfunc ( "stdH" , "sensH" , answer="YES" )
                    sensfunc ( "stdK" , "sensK" , answer="YES" )
                }
                if (noJ) {
                    imcombine ( procdir+outname+"_K_combine,"+procdir+outname+"_H_combine_m", procdir+outname+"_HK_combine" , combine="sum" , offset="wcs" )
                    if (background) {
                        background ( procdir+outname+"_HK_combine" , procdir+outname+"_HK_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background, sample=sample )
                    }
                }
                else {
                    imcombine ( procdir+outname+"_K_combine,"+procdir+outname+"_H_combine_m,"+procdir+outname+"_J_combine_m" , procdir+outname+"_JHK_combine" , combine="sum" , offset="wcs" )
                    if (background) {
                        background ( procdir+outname+"_JHK_combine" , procdir+outname+"_JHK_backsub" , axis=2 , naverage=2 , order=2 , high=2 , low=2 , niter=3 , interactive=interactive_background, sample=sample )
                    }
                }
            }
        }
    }

end

