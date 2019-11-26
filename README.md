HQUAILS
========
*Handy QUAsar emissIon Line fitS (pronounced Quails) by Raphael Hviding*
-------------

HQUAILS is a Python code designed to fit emission lines in the spectra of active galactic nuclei. In particular, it was built in order to fit AGN spectra where many of the parameters of the emission lines are tied with respect to one another. HQUAILS attempts to automate this process. For example, tying the redshifts and velocity widths of AGN lines (e.g. OIII, NII) together, but keeping that seperate from the redshifts and velocity widths of galaxy lines (e.g. Balmer series lines).

HQUAILS was also built in order to test the inclusion of additional fitting parameters. For example, is the spectrum better fit with a broad Halpha component? Or an outflowing OIII component? HQUAILS builds a base model based on the spectrum, and iteratively tests whethere different additional components are justified to add to the model, based on an F-test.

Requirements
-------------

HQUAILS was developed using Astroconda.
Instructions to install Astroconda can be found at:
https://astroconda.readthedocs.io/en/latest/

Express dependencies at time of writing are:
* Python 3.6.9
* NumPy 1.17.3
* AstroPy 3.2.3
* SciPy 1.3.1
* json 2.0.9

Installation
-------------

HQUAILS can be installed by cloning this Git repository. 

In your working directory, you need to copy "matplotlibrc" file to control the plotting settings. This is most important if you are running HQUAILS with multiprocessing as this file sets the matplotlib backed to "Agg", a non-interactive backend, required for generating and saving figures on multiple threads.

How it works
-------------

Parameter File
-------------

HQUAILS is controlled entirely by the "PARAMS.json" file. And example parameter file is included in the repository.

* Outfolder: This parameter is the path to the output directory. 
* RegionWidth: The border around emission lines (same units as spectrun wavelength).
* BackgroundDeg: Degree of polynomial for continuum background.
* MaxIter: Maximum number of minimazation algorithm iterations.
* NBoot: Number of bootstrap iterations to constrain error on parameters.
* FThresh: F-test threshold to incorproate additional model parameters.
* NPool: Number of processes to open with python multiprocessing. Set equal to 1 to use only a single thread.
* Plotting: Produce plots or not.
* EmissionLines: Dictionary of emission lines to be fit by HQUAILS. The structure of this dictionary is crucial to the operation of HQUAILS. The following section details the format of this dictionary.

Emission Line Dictionary
-------------

Here we describe the format of the emission line dictionary.

1. The emission line dictionary is made up of Groups. All spectral features in the same Group will share a common redshift. This means, during fitting, their redshifts will be forcibly tied to be equal.
2. Each Group is made out of Species. All spectral features in the same Species will share a velocity dispersion. This means, during fitting, their velocity dispersions will be forcibly tied to be equal.

   Each Species is also associated with an integer flag. The integer flag controls will additional parameters HQUAILS will attempt to add to the spectral features of this species. The value in each bit is a boolean flag for each kind of additional component, which can be found in the "AdditionalComponents.py" file. 
   
   Each additional component must be associated with a Group, which may or may not be the same Group as the original species. An additional list is passed to each species, specifying where each additional component that will attempt to be added must go. It can even be a group that does not exist in the emission line dictionary, it will simply be added. (Note: This means the sum of the binary integer flag must be equal to the length of the FlagGroups list.)


Each Species is assocaited with an integer flag. This flag decides which additional parameters HQUAILS will attempt to add to this species. The value of the bit in each position of the integer decides what type of extra component will attempt to be added. The  contains the types of additional components that can be added, and makes it easy to add more in the future. Finally, each species is associated with a list of Groups
3. Each species contains "lines". Lines are made up of a (rest) wavelength (same units as spectrun wavelength) and relative flux. Lines within each species will have their flux scaled relative to each other based on this factor.

The "PARAMS.json" file in the directory gives a good example of how to take advantage of these features. It consists of two groups:
1. AGN. These features will all share the same redshift. It is made up of two species. 
   1) OIII. These features will share the same velocity dispersion. In addition, these have been flagged
2. Galaxy

Running HQUAILS
-------------


HQUAILS can be run from any directory as long as the full path to the code is provided. 

HQUAILS ships with scripts for running:

1. "run_HQUAILS_single.py"

   This script is designed to run HQUAILS over a single object. 

HQUAILS
-------------
* CustomModels.py

  Here are where the custom models used in HQUAILS are defined. Here exists a gaussian emission line model and a polynomial continuum background. While the parameters are in the rest frame, the model fits in the observed frame.



License
-------------
HQUAILS is an open-source software available under the GNU General Public Licence. In a nutshell, this code can be used and distributed by anyone, but any code that includes HQUAILS must also be distributed freely and openly (see LICENCE file for details).

FAQ
-------------
**Why is it spelled HQUAILS but pronounced Quails?**

*The author of the code, R. E. Hviding (pronounced VEE-ding) thought it important to draw attention to names that start with an H where the H is not pronounced.*

**How can I fit my spectra with Voight or Lorentzian profiles?**

*In order to fit with with other functions, you can simply add them to the "CustomModels.py" file. However, you'll have to make sure they are built in the same framwork as the other models there.*