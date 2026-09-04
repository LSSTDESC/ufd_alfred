# Ultra Faint Dwarf Galaxies for Dark Matter

Pipeline tools for finding ultra-faint dwarf galaxies for dark matter analysis using joint LSST and Euclid photometry+morphology.  

:construction: under construction

setup:
- copy a venv with the lsst pipelines
    - e.g. `conda create --name alfred_venv --clone lsst-scipipe-12.3.0-exact`
    - `conda activate ufd_alfred`
- check butler and geom are there
    - if not, install lsst pipelines
    - maybe this works?:
        - `eups distrib install -t v30_0_10 lsst_distrib`
        - `setup lsst_distrib`
- `python -m pip install git+https://github.com/astropy/astroquery.git`
- `git clone https://github.com/DarkEnergySurvey/ugali.git && cd ugali`
- `python setup.py install`
- `git clone https://github.com/sidneymau/simple_adl.git`
    - so far just leaving the clone there 


Run by changing config.yaml to paths suitable for specific system and then run main.py.  
  
All pipeline functions live in alfred directory (named for the butler). Tests are to be written. Plots will be saved in the plots directory. Scratch is all my old code files. 
This currently uses fork/branch keexcell/simple\_adl/tree/euclid-codes for isochrones (uses an isochrone with Euclid and Roman bands available)  
 
Works on NERSC right now (some code in scratch was written to work on USDF. I believe the only thing that changes are the paths to data/collections).  
