# TspecCubes

This repository contains scripts used by Youngblood, Ginsburg, and Bally 2016 to create position-position-velocity cubes from APO TripleSpec slit scans and extract integrated intensity maps of emission lines and kinematic information from them.

Contents (in order of use):

-- pipeline_jhk.cl -- this IRAF script takes raw TripleSpec spectra, crops the JHK orders, spatially and spectrally rectifies the individual orders, and stitches them together. This is modified from pipeline_jhk.cl in https://github.com/keflavich/agpy.

-- reduction_script.txt -- this text file shows the user how prepare the spectra for pipeline_jhk.cl.

-- process_jhk.py -- Python script that flat-fields, airglow and background subtracts, flux calibrates, and applies a radial velocity correction.

-- make_cubes.py -- Python script that stitches the processed JHK spectra into a data cube with WCS coords.

-- combine_cubes.py -- Python script that combines the individual cubes into one big cube.

-- fit_cube.py -- Python script to create moment maps of an emission line from the big cube.
