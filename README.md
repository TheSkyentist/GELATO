HQUAILS
========
*Handy QUAsar emissIon Line fitS (pronounced Quails) by Raphael Hviding*
-------------

HQUAILS is a Python code designed to fit emission lines in the spectra of active galactic nuclei. In particular, it was built in order to fit AGN spectra where many of the parameters of the emission lines are tied with respect to one another. HQUAILS attempts to automate this process. For example, tying the redshifts of AGN lines (e.g. OIII, NII) together, and the flux ratios of the lines therein, but keeping that seperate from the redshifts of galaxy lines (e.g. Balmer series lines).

HQUAILS was also built in order to test the inclusion of additional fitting parameters. For example, is the spectrum better fit with a broad Halpha component? Or an outflowing OIII component? HQUAILS builds a base model based on the spectrum, and iteratively tests whethere different additional components are justified to add to the model, based on an F-test.

HQUAILS was designed to be run on SDSS spectra using a LM non-linear least squares algorithm with Gaussian line profiles. However, it is written in such a way that these assumptions can be switched out as needed. Please read through the entire documentation to see how this can be done.

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

1. First, the spectrum is loaded. Here, based on the emission line dictionary and redshift provided, the code determines which emission lines actually lie inside the domain of the spectrum. It then constructs regions around these emission lines based on the region width provided. If regions overlap, they are merged.

2. The base model is then constructed based on the emission line dictionary. The starting values are generated based on the spectrum. The model is then fit to the spectrum. The default fitting minization is the Levenberg–Marquardt non-linear least squares algorithm. This can be adjusted.

3. The additional components are then added to the base model and tested seperately. If the fit is statistically better with the additional component, it is accepted. This is decided by performing and F-test. All accepted additional components are then collected and incorporated into the final model.

4. In order to constraint fit uncertainites, the flux is bootstrapped with respect to provided uncertainties and the fit is run again. This process is repeated as many times as required by the user. 

5. The full set of bootstrapped parameters is then saved to disk. Finally, a figure of the final fit is produced and saved. There exists a convinience function for finding the median values of each spectrum model fit and collecting them into one final table.

Parameter File
-------------

The behaviour of HQUAILS is controlled entirely by the "PARAMS.json" file. And example parameter file is included in the repository.

* Outfolder: This parameter is the path to the output directory. 
* RegionWidth: The border around emission lines (same units as spectrum wavelength).
* BackgroundDeg: Degree of polynomial for continuum background.
* MaxIter: Maximum number of minimazation algorithm iterations.
* NBoot: Number of bootstrap iterations to constrain error on parameters.
* FThresh: F-test threshold to incorproate additional model parameters.
* NProcess: Number of processes to open with python multiprocessing. Set equal to 1 to use only a single thread.
* Plotting: Produce plots or not.
* Concatenate: To concatenate the results of a multiple HQUAILS run or not.
* EmissionLines: Dictionary of emission lines to be fit by HQUAILS. The structure of this dictionary is crucial to the operation of HQUAILS. The following section details the format of this dictionary.

Emission Line Dictionary
-------------

Here we describe the format of the emission line dictionary.

1. The emission line dictionary is made up of Groups. All spectral features in the same Group will share a common redshift. This means, during fitting, their redshifts will be forcibly tied to be equal.

2. Each Group is made out of Species. All spectral features in the same Species will share a velocity dispersion. This means, during fitting, their velocity dispersions will be forcibly tied to be equal.

   Each Species is also associated with an integer flag. The integer flag controls will additional parameters HQUAILS will attempt to add to the spectral features of this species. The value in each bit, from right to left, is a boolean flag for each kind of additional component, which can be found in the "AdditionalComponents.py" file. 
   
   Each additional component must be associated with a Group, which may or may not be the same Group as the original species. An additional list is passed to each species, specifying where each additional component that will attempt to be added must go. It can even be a group that does not exist in the emission line dictionary, it will simply be added. (Note: This means the sum of the binary integer flag must be equal to the length of the FlagGroups list.)

3. Each Species is made out of Lines. Lines are made up of a (rest) wavelength (same units as spectrum wavelength) and relative flux. Lines within each species will have their flux scaled relative to each other based on this factor. This means, during fitting, their fluxes will be tied to have the given relative values.

The "PARAMS.json" file in the directory gives a good example of how to take advantage of these features. It consists of two groups:

1. AGN. These features will all share the same redshift. It is made up of two species. 
   1. OIII. These features will share the same velocity dispersion. In addition, these have been flagged with a 1, which corresponds to an larger velocity component. This additional component will be placed in its own new group, which will be labelled "Outflow". It is made out of two lines.
      * A line with a rest wavelength of 5006.77 and a relative flux of 1. 
      * A line with a rest wavelength of 4958.83 and a relative flux of 0.35. This means this line will always have 0.35/1 times the flux of the first line.
   2. NII. These features will share the same velocity dispersion. There are no flags on this component, so the list is empty. It is made out of two lines.
      * A line with a rest wavelength of 6583.34 and a relative flux of 1. 
      * A line with a rest wavelength of 6547.96 and a relative flux of 0.34. This means this line will always have 0.34/1 times the flux of the first line.
2. Galaxy
   1. Halpha. These have been flagged with a 1, which corresponds to an larger velocity component. This additional component will be placed in its own new group, which will be labelled "Broad". It is made out of one line.
      * A line with a rest wavelength of 6562.80 and a relative flux of 1. Since there is only one line, the relative flux value does not matter.
   2. Hbeta. These have been flagged with a 3, or in binary, 11. This corresponds to both a larger velocity component and an absorption component. These will be placed in the "Broad" and "Absorption" groups respetively. Not that if both "Broad" lines are accepted (from Hbeta and Halpha), they will share the same redshift by design. It is made out of one line.
      * A line with a rest wavelength of 4861.32 and a relative flux of 1. Since there is only one line, the relative flux value does not matter.

Running HQUAILS
-------------

In order to run HQUAILS you need:

* The PARAMS.json file. 
* The spectrum or spectra.
* The redshift of each spectrum. The redshift of the object must be passed to construct the spectrum object. While the redshift is a fitted parameter, the provided value must be correct to at least 1 part in 100. A basic estimate from the apparent position of any identified emission line should suffice.
* (If plotting) the matplotlibrc file in your working directory, especially if you are running on multiple threads, in which case the non-interactive backend must be specified. 

The two wrappers for HQUAILS are:

1. "run_HQUAILS_single.py"

   This script is designed to run HQUAILS over a single object. This takes 3 positional arguements, the path to the parameters file, the path to the spectrum, and the redshift of the object. 

```
python ~/Documents/HQUAILS/run_HQUAILS_multi.py ~/Example/PARAMS.json ~/Data/spectrum.fits 1.122
```

2. "run_HQUAILS_multi.py"

   This script is designed to run HQUAILS over a list of objects. This takes 2 positional arguements, the path to the parameters file, and the path to the list of objects. 

```
python ~/Documents/HQUAILS/run_HQUAILS_multi.py ~/Example/PARAMS.json ~/Data/spectra_with_redshifts.txt
```

The plots for HQUAILS can also be created directly from the spectra and the results file in the following manners:

For a single plot:

```
python ~/Documents/HQUAILS/Plotting.py ~/Example/PARAMS.json --Spectrum ~/Data/spectrum.fits --Redshift 1.122
```

For multiple plots:

```
python ~/Documents/HQUAILS/run_HQUAILS_multi.py ~/Example/PARAMS.json --ObjectList ~/Data/spectra_with_redshifts.txt
```

The concated results for HQUAILS can also be created directly the results files in the following manners:

```
python ~/Documents/HQUAILS/run_HQUAILS_multi.py ~/Example/PARAMS.json ~/Data/spectra_with_redshifts.txt
```

HQUAILS cast (in order of appearance)
------------
* README.md
  
  Here you are! The documenation for HQUAILS.

* run_HQUAILS_single.py
  
  Wrapper for running HQUAILS on a single object. 
  
* run_HQUAILS_multi.py
  
  Wrapper for running HQUAILS on multiple objects. If specifying mulitple processes, each object will be run on an independant thread. To load an object file differently, this file should be edited.

* PARAMS.json
  
  Paramter file that controls HQUAILS behaviour.

* ConstructParams.py

  Routines for turning the PARAMS.json file into a python dictionary, and verifying that it is in the correct format. 

* HQUAILS.py

  Main HQUAILS function that calls and coordinates the whole operation.

* SpectrumClass.py

  A class the defines how a spectrum is loaded. HQUAILS was designed for SDSS spectra. To load in any other kind of spectrum, you can edit the way this class is initalized.

* BuildModel.py

  Thise file handles the construction of models and for tying the various parameters together as outlined by the emission line dictionary.

* CustomModels.py

  Here are where the custom models used in HQUAILS are defined. Here exists a gaussian emission line model and a polynomial continuum background. The parameters for each model are defined with respect to the rest frame, but the output of the model is in the observed frame.

* AdditionalComponents.py

  Here are where the additional components are defined along with their bit flag positions. In order to add extra additional components, this file can be easily extended to include more models.

* FittingModel.py

  Here are the scripts for fitting HQUAILS generated models and for testing the inclusion of additional parameters. To change the fitting algorithm, this file can be edited. 

* FischerTest.py

  Here are scripts for performing an F-test to test whethere the inclusion of an additional model parameter is statistically better.

* Plotting.py

  Here are the scripts for creating and saving figures of the fits. Can also be run directly on HQUAILS results in order to create figures after the fact. 

* ConcatenateResults.py

  Scripts for contatenating results from a multi HQUAILS run. Can also be run independantly on results after the fact. 

* matplotlibrc

  A matplotlib settings file that controls the output of figures. Can be changed to the user's liking. However, for running HQUAILS on multiple threads, the backend must be set to a non-interactive backend, e.g. "Agg". 

* LICENSE

  Code license, HQUAILS is distributed under the GNU General Public Licence 3.

License
-------------
HQUAILS is an open-source software available under the GNU General Public Licence 3. In a nutshell, this code can be used and distributed by anyone, but any code that includes HQUAILS must also be distributed freely and openly (see LICENCE file for details).

FAQ
-------------
**Why is it spelled HQUAILS but pronounced Quails?**

*The author of the code, R. E. Hviding (pronounced VEE-ding) thought it important to draw attention to names that start with an H where the H is not pronounced.*

**How can I fit my spectra with Voight or Lorentzian profiles?**

*In order to fit with with other functions, you can simply add them to the "CustomModels.py" file. However, you'll have to make sure they are built in the same framwork as the other models there.*