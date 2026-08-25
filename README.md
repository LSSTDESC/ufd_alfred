# Ultra Faint Dwarf Galaxies for Dark Matter

Pipeline tools for finding ultra-faint dwarf galaxies for dark matter analysis using joint LSST and Euclid photometry+morphology.  
Run by changing config.yaml to paths suitable for specific system and then run main.py.  
  
All pipeline functions live in alfred directory (named for the butler). Tests are to be written. Plots will be saved in the plots directory. Scratch is all my old code files. 
This currently uses fork/branch keexcell/simple\_adl/tree/euclid-codes for isochrones (uses an isochrone with Euclid and Roman bands available)  
 
Works on NERSC right now (some code in scratch was written to work on USDF. I believe the only thing that changes are the paths to data/collections).  
