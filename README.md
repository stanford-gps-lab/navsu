# navsu

NavSU is a MATLAB toolbox for general GNSS processing.  Currently contains parsers, time functions, tools to download and and handle IGS products, and more.  

# How to run:
Just put the +navsu toolbox in your path!  See the \examples folder for run examples.  

# Sub-toolboxes:
## svOrbitClock
Class for downloading and handling GNSS corrections, including orbit, clock, ionospheric, antenna phase center, and differential code biases.  See the example \examples\example_setup_svOrbitClock.m for some more information.

## ftp
Main function here: navsu.ftp.download to download various IGS products.

## readfiles
Lots of parsing tools here.  They can be called directly from here, or they can be used by the svOrbitClock object. 
To parse RINEX obsevation files: navsu.readfiles.loadRinexObs
To parse RINEX orbit .sp3 files: navsu.readfiles.readSp3

## svprn
Tools for converting to and from SVN and PRN as well as additional information about satellites including block and frequency assignements (for GLONASS)

## time
Time utilities, convert to and from GPS epochs (seconds since start of GPS time), julian date, MATLAB datenum, and calendar date. 



