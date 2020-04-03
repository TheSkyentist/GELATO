HQUAILS
========
*Handy QUAsar emissIon Line fitS (pronounced Quails) by Raphael Hviding*
-------------

HQUAILS is a Python code designed to fit emission lines in the spectra of active galactic nuclei. In particular, it was built in order to fit AGN spectra where many of the parameters of the emission lines are tied with respect to one another. HQUAILS attempts to automate this process. For example, tying the redshifts of AGN lines (e.g. OIII, NII) together, and the flux ratios of the lines therein, but keeping that separate from the redshifts of galaxy lines (e.g. Balmer series lines).

HQUAILS was also built in order to test the inclusion of additional fitting parameters. For example, is the spectrum better fit with a broad Halpha component? Or an outflowing OIII component? HQUAILS builds a base model based on the spectrum, and iteratively tests whether different additional components are justified to add to the model, based on an F-test.

HQUAILS was designed to be run on SDSS spectra using a LM non-linear least squares algorithm with Gaussian line profiles. However, it is written in such a way that these assumptions can be switched out as needed. Please read through the entire documentation to see how this can be done.

Requirements
-------------

HQUAILS was developed using Astropy 3.2.3 and Python 3.6.10.

To install the dependancies, I recommend installing conda (through [Miniconda](https://docs.conda.io/en/latest/miniconda.html)).

The environment can be installed from the provided "environment.yml" file.

```bash
conda env create -f environment.yml
```

The environment will then be installed under the name HQUAILS and can then be activated.

```bash
conda activate HQUAILS
```

Whenever running HQUAILS scripts, it must be run from this environment.

Installation
-------------

HQUAILS can be installed by cloning this git repository.

In your working directory, **you need to copy the "matplotlibrc" file** to control the plotting settings. This is most important if you are running HQUAILS with multiprocessing as this file sets the matplotlib backed to "Agg", a non-interactive backend, required for generating and saving figures on multiple threads.

How it works
-------------

1. First, the spectrum is loaded. Here, based on the emission group dictionary and redshift provided, the code determines which emission lines actually lie inside the domain of the spectrum. It then constructs regions around these emission lines based on the region width provided. If regions overlap, they are merged.

2. The base model is then constructed based on the emission line dictionary. The starting values are generated based on the spectrum. The model is then fit to the spectrum. The default fitting minimization is the Levenberg–Marquardt non-linear least squares algorithm. This can be adjusted.

3. The additional components are then added to the base model and tested separately. If the fit is statistically better with the additional component, it is accepted. This is decided by performing and F-test. The combinations of all accepted additional components are then then tested by measuring their AICs. The model set with the lowest AIC is the final model.

4. In order to constraint fit uncertainties, the flux is bootstrapped with respect to provided uncertainties and the fit is run again. This process is repeated as many times as required by the user.

5. The full set of bootstrapped parameters is then saved to disk. Finally, a figure of the final fit is produced and saved. There exists a convenience function for finding the median values of each spectrum model fit and collecting them into one final table.

6. The median of each parameter is then found for each object, and the results from each object are then joined into one final results file.

Parameter File
-------------

The behaviour of HQUAILS is controlled entirely by the "PARAMS.json" file. And example parameter file is included in the repository.

* Outfolder: This parameter is the path to the output directory. 
* RegionWidth: The border around emission lines in velocity space that is fit (km/s).
* LineDataWidth: The border around an emission line that must be contained within the spectrum in order to be fit, also where the flux is calculate to assign an initial value.
* BackgroundDeg: Degree of polynomial for continuum background.
* MaxIter: Maximum number of minimization algorithm iterations.
* NBoot: Number of bootstrap iterations to constrain error on parameters.
* FThresh: F-test threshold to incorporate additional model parameters.
* NProcess: Number of processes to open with python multiprocessing. Set equal to 1 to use only a single thread.
* Plotting: Produce plots or not.
* PlotComp: To plot only the components of the fit or the total fit.
* Concatenate: To concatenate the results of a multiple HQUAILS run or not.
* Verbose: To print output for each object or not.
* EmissionGroups: Dictionary of emission lines to be fit by HQUAILS. The structure of this dictionary is crucial to the operation of HQUAILS. The following section details the format of this dictionary.

Emission Line Dictionary
-------------

Here we describe the format of the emission line dictionary.

1. The emission groups dictionary is made up of Groups. All spectral features in the same group can share a common redshift, a common dispersion, or neither.  This means, during fitting, their redshifts or dispersions can be forcibly tied to be equal.

      * Each group has a Name, which controls how its parameters appear in the output.
      * Each group has TieRedshift flag, which controls if the redshifts of all the group elements are tied or not.
      * Each group has TieSigma flag, which controls if the dispersions of all the group elements are tied or not.
      * Finally, each group is made out of a list of species.

2. Each Group contains a list of Species. All spectral features in the same Species will share a redshift velocity and dispersion. This means, during fitting, their velocity dispersions and redshifts will be forcibly tied to be equal.

      * Each species has a name, which controls how its parameters appear in the output.
      * Each species has a Flag. The integer flag controls will additional parameters HQUAILS will attempt to add to the spectral features of this species. The value in each bit, from right to left (increasing order of magnitude), is a boolean flag for each kind of additional component, which can be found in the "AdditionalComponents.py" file. 
      * Each species has a FlagGroup. Each additional component must be associated with a Group, which may or may not be the same Group as the original species. An additional list is passed to each species, specifying where each additional component that will attempt to be added must go. The group must exist, even if empty, as it needs to be created with the flags. (Note: This means the sum of the binary integer flag must be equal to the length of the FlagGroups list.)

3. Each Species contains a list of lines. Lines can be set to have relative fluxes.  This means, during fitting, their fluxes will be tied to have the given relative values.

      * Lines have a Wavelength. This is the rest wavelength of the line (same units as spectrum wavelength).
      * Lines have a RelStrength. This is a relative strength to the other members of the species. If set to null, it will have an independent flux.

The "PARAMS.json" file in the directory gives a good example of how to take advantage of these features. It consists of five groups:

1. Name: NLR. Here these features have not been set to share redshifts nor dispersions. It has a list of species:
   1. Name: SII. These features will share the same velocity dispersion and redshift. There are no flags on this component, so the list is empty. It is made out of two lines.
      * A line with a rest wavelength of 6716.31 and a relative flux of null.
      * A line with a rest wavelength of 6730.68 and a relative flux of null. This means the line fluxes are completely independent.
   2. Name: NII. These features will share the same velocity dispersion and redshift. There are no flags on this component, so the list is empty. It is made out of two lines.
      * A line with a rest wavelength of 6583.34 and a relative flux of 1.
      * A line with a rest wavelength of 6547.96 and a relative flux of 0.34. This means this line will always have 0.34/1 times the flux of the first line.
   3. OIII. These features will share the same velocity dispersion and redshift. These have been flagged with a 2, or in binary 10. This corresponds to an outflowing velocity component. This additional component will be placed in its own new group, which will be labelled "Outflow". It is made out of two lines.
      * A line with a rest wavelength of 5006.77 and a relative flux of 1.
      * A line with a rest wavelength of 4958.83 and a relative flux of 0.35. This means this line will always have 0.35/1 times the flux of the first line.
   4. Name: [NeIII]. Singlet line with no flags.
   5. Name: [OII]. Singlet line with no flags.
   6. Name: [NeV]. Singlet line with no flags.
2. Name: Stars. Here these features are set to share redshifts and dispersions. It has a list of species:
   1. Halpha. A singlet line flagged with a 1, which corresponds to an larger velocity component. This additional component will be placed the group named "Broad". It is made out of one line.
      * A line with a rest wavelength of 6562.80 and a relative flux of null.
   2. Hbeta. These have been flagged with a 5, or in binary, 101. This corresponds to both a larger velocity component and an absorption component. These will be placed in the "Broad" and "Absorption" groups respectively. It is made out of one line.
      * A line with a rest wavelength of 4861.32 and a relative flux of 1. Since there is only one line, the relative flux value does not matter.
3. Name: Outflow. This is an empty group as it may receive additional components from other groups. If more than one component lands in this group, they will not share redshifts and dispersions.
4. Name: Broad. This is an empty group as it may receive additional components from other groups. If more than one component lands in this group, they will share redshifts and dispersions. E.g. if both "Broad" lines are accepted (from Hbeta and Halpha), they will share the same redshift and dispersion by design.
5. Name: Absorption. This is an empty group as it may receive additional components from other groups. If more than one component lands in this group, they will not share redshifts and dispersions.

Running HQUAILS
-------------

In order to run HQUAILS you need:

* The PARAMS.json file. 
* The spectrum or spectra.
* The redshift of each spectrum. The redshift of the object must be passed to construct the spectrum object. While the redshift is a fitted parameter, the provided value must be correct to at least 1 part in 100. A basic estimate from the apparent position of any identified emission line should suffice.
* (If plotting) the matplotlibrc file in your working directory, especially if you are running on multiple threads, in which case the non-interactive backend must be specified. 

The two wrappers for HQUAILS are:

1. "run_HQUAILS_single.py"

   This script is designed to run HQUAILS over a single object. This takes 3 positional arguments, the path to the parameters file, the path to the spectrum, and the redshift of the object.

  ```bash
  python ~/Documents/HQUAILS/run_HQUAILS_multi.py ~/Example/PARAMS.json ~/Data/spectrum.fits 1.122
  ```

2. "run_HQUAILS_multi.py"

   This script is designed to run HQUAILS over a list of objects. This takes 2 positional arguments, the path to the parameters file, and the path to the list of objects. 

```bash
python ~/Documents/HQUAILS/run_HQUAILS_multi.py ~/Example/PARAMS.json ~/Data/spectra_with_redshifts.txt
```

The plots for HQUAILS can also be created directly from the spectra and the results file in the following manners:

For a single plot:

```bash
python ~/Documents/HQUAILS/Plotting.py ~/Example/PARAMS.json --Spectrum ~/Data/spectrum.fits --Redshift 1.122
```

For multiple plots:

```bash
python ~/Documents/HQUAILS/run_HQUAILS_multi.py ~/Example/PARAMS.json --ObjectList ~/Data/spectra_with_redshifts.txt
```

The concatenated results for HQUAILS can also be created directly the results files in the following manners:

```bash
python ~/Documents/HQUAILS/run_HQUAILS_multi.py ~/Example/PARAMS.json ~/Data/spectra_with_redshifts.txt
```

Running the Example
-------------
Here are the following instructions to run HQUAILS. This tutorial assumes you start in the Example directory. First we need to activate our astroconda environment.

```bash
conda activate HQUAILS
```

We can then run the code over the whole data set.

```bash
python ../run_HQUAILS_multi.py ExPARAMS.json ExObjList.csv
```

In order to produce a concatenated table of all of the results, we can run the following code. (We could have achieved the same result by changing the Concatenate parameter to true)

```bash
python ../ConcatResults.py ExPARAMS.json ExObjList.csv
```

This will produce result tables and some plots in the Results folder. We can then edit the value of the PlotComp parameter and change it to false in order to produce the rest of the plots.

```bash
python ../Plotting.py ExPARAMS.json --ObjectList ExObjList.csv
```

HQUAILS cast (in order of appearance)
------------
* README.md
  
  Here you are! The documentation for HQUAILS.

* run_HQUAILS_single.py
  
  Wrapper for running HQUAILS on a single object. 
  
* run_HQUAILS_multi.py
  
  Wrapper for running HQUAILS on multiple objects. If specifying multiple processes, each object will be run on an independent thread. To load an object file differently, this file should be edited.

* PARAMS.json
  
  Paramter file that controls HQUAILS behavior.

* ConstructParams.py

  Routines for turning the PARAMS.json file into a python dictionary, and verifying that it is in the correct format. 

* HQUAILS.py

  Main HQUAILS function that calls and coordinates the whole operation.

* SpectrumClass.py

  A class the defines how a spectrum is loaded. HQUAILS was designed for SDSS spectra. To load in any other kind of spectrum, you can edit the way this class is initialized.

* BuildModel.py

  This file handles the construction of models and for tying the various parameters together as outlined by the emission line dictionary.

* CustomModels.py

  Here are where the custom models used in HQUAILS are defined. Here exists a gaussian emission line model and a polynomial continuum background. The parameters for each model are defined with respect to the rest frame, but the output of the model is in the observed frame. This is where the velocity width limits on emission features can bs set. 

* AdditionalComponents.py

  Here are where the additional components are defined along with their bit flag positions. In order to add extra additional components, this file can be easily extended to include more models. This is where the velocity dispersion limits on additional components can be modified.

* FittingModel.py

  Here are the scripts for fitting HQUAILS generated models and for testing the inclusion of additional parameters. To change the fitting algorithm, this file can be edited.

* ModelComparison.py

  Here are scripts for model comparison and selection, including F-Tests and AIC calculation.

* Plotting.py

  Here are the scripts for creating and saving figures of the fits. Can also be run directly on HQUAILS results in order to create figures after the fact. The plots can either be plotted as the components or the final result only.

* ConcatenateResults.py

  Scripts for concatenating results from a multi HQUAILS run. Can also be run independently on results after the fact.

* matplotlibrc

  A matplotlib settings file that controls the output of figures. Can be changed to the user's liking. However, for running HQUAILS on multiple threads, the backend must be set to a non-interactive backend, e.g. "Agg".

* LICENSE

  Code license, HQUAILS is distributed under the GNU General Public License 3.

License
-------------
HQUAILS is an open-source software available under the GNU General Public License 3. In a nutshell, this code can be used and distributed by anyone, but any code that includes HQUAILS must also be distributed freely and openly (see LICENSE file for details).

FAQ
-------------
**Why is it spelled HQUAILS but pronounced Quails?**

*The author of the code, R. E. Hviding (pronounced VEE-ding) thought it important to draw attention to names that start with an H where the H is not pronounced.*

**How can I fit my spectra with Voight or Lorentzian profiles?**

*In order to fit with with other functions, you can simply add them to the "CustomModels.py" file. However, you'll have to make sure they are built in the same framework as the other models there.*

**What are the units?**

*As long as your RegionWidth, line Wavelengths, and spectrum wavelengths are all in the same units, the fitting will work. Dispersions are given in km/s. The units in plotting can be changed in the Plotting.py file.*

**Do you mean velocity offsets, not redshifts?**

*Each emission line is characterized by a redshift, which is trivial to convert to a velocity offset once a reference line is chosen. However this requires the user to choose a reference line. HQUAILS remains agnostic to this procedure and simply returns the redshift of each line.*
