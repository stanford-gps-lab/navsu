# navsu

NavSU is a MATLAB toolbox for general GNSS processing.  Currently contains parsers, time functions, tools to download and and handle IGS products, and more.  

# How to run:
Just put the +navsu toolbox in your path!  See the \examples folder for run examples.  

# Sub-toolboxes:
## svOrbitClock
Class for downloading and handling GNSS corrections, including orbit, clock, ionospheric, antenna phase center, and differential code biases.  See the example \examples\example_setup_svOrbitClock.m for some more information.

## ftp
Main function here: navsu.ftp.download to download various IGS products.

## lsNav
Contains a conventional least squares navigation engine. Capable of providing dual frequency, multi constellation solutions. Designed to offer flexibility in used frequencies, constellations and signals. Carrier smoothing is accomplished by separate class. See examples/example_lsNavEngine.m.

## readfiles
Lots of parsing tools here.  They can be called directly from here, or they can be used by the svOrbitClock object. 
To parse RINEX obsevation files: navsu.readfiles.loadRinexObs
To parse RINEX orbit .sp3 files: navsu.readfiles.readSp3

## svprn
Tools for converting to and from SVN and PRN as well as additional information about satellites including block and frequency assignements (for GLONASS)

## time
Time utilities, convert to and from GPS epochs (seconds since start of GPS time), julian date, MATLAB datenum, and calendar date. 

## Matlab Continuous Integration

This build checker does the following.
1. Executes Matlab's `checkcode` function on all repository files and asserts no suggestions.
2. Finds all Matlab tests in the repository, executes Matlab's `runtests`, then reports the results.
3. If their is a pull request associated with the push, a unit test coverage report is attached to the report.

This build checker will also perform standard checks on repository Matlab code using Python with MISS_HIT.
Source can be found [here](https://github.com/florianschanda/miss_hit).
Documentation and installation instructions can be found [here](https://florianschanda.github.io/miss_hit).

To run the all the checks locally, one must execute the following from bash at the repository root directory (note the `--fix` will automatically fix style issues)
```bash
mh_style --process-slx --fix
mh_metric --ci
mh_lint
```
and the following from Matlab at the repository root directory (with the appropriate Matlab path setup).
```matlab
addpath('.github/workflows');
matlab_checkcode_on_directory('.');
matlab_runtests_on_directory('.');
```




