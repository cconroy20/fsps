FSPS: Flexible Stellar Population Synthesis
=====
version 3.2

References
---------
When using this code please cite the following papers:
 * Conroy, Gunn, & White 2009, ApJ, 699, 486
 * Conroy & Gunn 2010, ApJ, 712, 833

Installation
----------
If you have git installed, FSPS can be obtained with the following commands:
```
cd /path/to/desired/location/
git clone https://github.com/cconroy20/fsps
```
Otherwise download a gzipped tarball from [here](https://github.com/cconroy20/fsps/releases). Then follow the instructions at [doc/INSTALL](doc/INSTALL).

You should not need to update the git repository until an update is announced (which is why you need to be on the mailing list - see [doc/INSTALL](doc/INSTALL)).  If you've obtained FSPS using git then when an update is announced you will need to simply type ``cd $SPS_HOME; git pull`` and then recompile.  If you have made your own edits to the FSPS files, git will attempt to gracefully merge your local version with the repository version.

Documentation
------
See the [Manual](doc/MANUAL.pdf)

Contents
---------
Below is a brief description of the contents of the directories in the
fsps root directory:

 * `ISOCHRONES`: Contains the isochrone tables for the BaSTI and Padova
isochrone sets.  The Geneva isochrones have been pasted onto the BaSTI
and Padova tables for high masses (M>70Msun), and the low-mass Lyon
models have been pasted on at low masses.  You should not edit these
files unless you know what you're doing.

 * `OUTPUTS`: Contains the outputs of a few example calls of the routines
autosps and simple.  You may wish to use this directory for all
outputs of the fsps routines.

 * `SPECTRA`: Contains the spectral libraries, the spectrum of an A0V star
used to set the Vega magnitude zero points, and a spectrum of the Sun.
The BaSeL spectra (based on the Kurucz models) are in binary format,
primarily to make the read in time faster and to decrease the size of
the fsps download.  The Hot_spectra directory contain the libraries
for O stars, WR stars, and post-AGB stars, from Smith et al. 2002 and
Rauch 2003, respectively.

 * `data`: Contains files that define the set of filters and indices used
in FSPS and the tabulated imfs and sfhs if those options are set.  The
files in this directory are readily user editable.

 * `doc`: Contains the manual, revision history, and installation
instructions.

* `dust`: Contains the dust attenuation curves for the Witt & Gordon
(2000) dust model and the dust emission spectra from the Draine & Li
2007 grain model.  Also contains the circumstellar dust models from 
Villaume et al. 2015 and the AGN dusty torus models of Nenkova et al. 2008.

* `nebular`: Contains the Cloudy lookup tables for nebular emission 
(both continuum and line emission) computed by Nell Byler.

* `pro`: Contains IDL files for reading in the .mag, .indx, and .spec
output files

* `src`: Contains the source files and routines from Numerical Recipes.






