GELATO
========

![Logo](./Images/Logo.png)
![Example of Fit](./Images/Example1.png)
![Example of Components](./Images/Example2.png)

*Galaxy/AGN Emission Line Analysis TOol by Raphael Hviding*
-------------

GELATO is a Python code designed to fit emission lines in the spectra of star forming galaxies and active galactic nuclei. In particular, it was built in order to fit spectra where many of the parameters of the emission lines are tied with respect to one another. GELATO attempts to automate this process. For example, tying the redshifts of AGN lines (e.g. OIII, NII) together, and the flux ratios of the lines therein, but keeping that separate from the redshifts of galaxy lines (e.g. Balmer series lines). In addition, GELATO is designed to fit additional components to lines with comples kinematics. For example, is the spectrum better fit with a broad Halpha component? Or an outflowing OIII component? GELATO builds a base model based on the spectrum, and iteratively tests whether different additional components are justified to add to the model, based on an [F-test](https://en.wikipedia.org/wiki/F-test) and then comparisons of [Akaike Information Criteria](https://en.wikipedia.org/wiki/Akaike_information_criterion).

The spectra are fit using a [Trust Region Reflective](https://epubs.siam.org/doi/10.1137/S1064827595289108) bounded non-linear least-squares optimization algorithm with Gaussian line profiles.


Installation
-------------

First, clone the GELATO git repository. If you have not used git before, this is easily done by using the commands before and cloning the directory over HTTPS.

```bash
cd /path/to/insallation/directory
git clone https://github.com/TheSkyentist/GELATO.git
```

GELATO is built primarily using Python. It primarily uses NumPy for math, SciPy for optimization, Astropy for FITS handling, matplotlib for plotting, and Jupyter for the example notebook.

To install the dependancies, I recommend installing conda (through [Miniconda](https://docs.conda.io/en/latest/miniconda.html)).

A conda environment with early all the dependencies can be installed via the provided "environment.yml" file. If you do not wish to use conda, you can install the dependencies enumerated in the "environment.yml" file.

```bash
cd /path/to/GELATO/directory
conda env create -f environment.yml
```

The environment will then be installed under the name GELATO and can then be activated.

```bash
conda activate gelato
```

Whenever running GELATO scripts, they must be run from this environment.

Then the GELATO scripts can be installed. Make sure you are in the GELATO conda environment.

```bash
cd /path/to/GELATO/directory
conda activate gelato
python setup.py install
```

In your working directory, **you need to copy the "matplotlibrc" file** to control the plotting settings. This is most important if you are running GELATO with multiprocessing as this file sets the matplotlib backed to "Agg", a non-interactive backend, required for generating and saving figures on multiple threads.

Updating GELATO
-------------

In order to update GELATO, you need to update the core packages, pull from the repo, and reinstall GELATO to your path. 

```bash
conda update -n gelato --all
cd /path/to/GELATO/directory
git pull
conda activate gelato
python setup.py install
```

In the case where a new version of GELATO releases where there are more strict dependecies, you may need to fully delete GELATO and its associated conda environment and reinstall from scratch.

How it works
-------------

1. Gathering Ingredients: First, the spectrum is loaded. The code assumes the spectrum file is a FITS table with the following columns and column names:
    1. The log10 of the wavelengths in Angstroms, column name: "loglam"
    2. The spectral flux density in flam units, column name: "flux"
    3. The inverse variances of the data points, column name: "ivar"

    Based on the emission line dictionary and redshift provided, the code determines which emission lines actually lie inside the domain of the spectrum. The region free from emission lines is then determined which will be used to obtain the initial fit to the continuum.

2. Creating Base (Continuum): GELATO models the continuum as a combination of Simple Stellar Populations (SSPs) from the [Extended MILES stellar library](http://research.iac.es/proyecto/miles/). We take SSP models assuming a Chabrier IMF (slope=1.3), the isochrones of Girardi et al. (2000) (Padova+00) with solar alpha abundance, and spanning a range of representatives metallicities and ages ([M/H] = [-1.31, -0.40, 0.00] and Age = [00.0631, 00.2512, 01.0000, 04.4668, 12.5893] (Gyr)) with nominal resolutions of 5 Angstroms. Note, since the continuum models have a minimum wavelength of 1680 Angstroms, there is a maximum wavelength that can be fit with GELATO based on the spectral coverage of the input spectra. The redshift is allowed to vary around the input redshift and the SSP models are fit to the region of continuum free from emission lines. The coefficients for the SSP models are constrained to be positive. Following the initial fit, an additional power law component is added, required to have a negative power law index and a positive coefficient. If the continuum model with a power law passes an F-test for its inclusion, it is added to the model. Finally, the redshift of the continuum model is frozen and not moving forward.

3. Creating Base (Emission Lines): The emission line models are then constructed based on the emission line dictionary. The starting values are generated based on the spectrum by looking at the range of values where the emission line would be expected to lie. The model flux is reasonable bounded based on these values, and the redshift of the line is bounded to be within 0.005 of it's starting value. The model is then fit to the spectrum.

4. Adding Flavor: The additional components are then added to the base model and tested separately. If the fit is statistically better with the additional component, it is accepted. This is decided by performing an F-test. The combinations of all accepted additional components are then then tested by measuring their Akaike Information Criteria (AICs). The model set with the lowest AIC is the final model.

5. Scooping Portions: In order to constraint fit uncertainties, the flux is bootstrapped with respect to provided uncertainties and the fit is run again.

6. Presenting gelato: Figures depicting the fit, for the entire spectrum and zoomed into specific lines, are then saved to disk.

7. Measuring texture: From the results, the rest equivalent width for each emission line is calculated. The height of the continuum is found by taking the median continuum in a region around the emission line.

8. Freezing results: The full set of bootstrapped parameters are saved to disk.

9. Combining gelato: If running on multiple objects, the median parameters and standard deviations for all of the fits are concatenated into one file and saved to disk.

Running GELATO
-------------

In order to run GELATO you need:

* The parameters file, a JSON file, the format of which is described below. An example JSON file is also provided.
* The spectrum or spectra. The log10 of the wavelength in Angstroms of the spectrum must be provided along with the flux in Flam units. The inverse variance of the fluxes, in corresponding units, must also be provided.
* The redshift of each spectrum. The redshift of the object must be passed to construct the spectrum object. While the redshift is a fitted parameter, the provided value must be correct to at least 1 part in 200, preferable 1 part in 1000. A basic estimate from the apparent position of any identified emission line should suffice.
* If running on a list of spectra, GELATO takes either...
  * a comma delimited file (ending in .csv) where each object occupies a different line. The first item in each line is the path to the spectrum. The second is the redshift of the spectrum.
  * a FITS table (ending in .fits) where each object occupies a different entry in the table. The table must have the column "Path" for the path to the object, and "z"
  containing the redshift of the object.
* (If plotting) the matplotlibrc file in your working directory, especially if you are running on multiple threads, in which case the non-interactive backend must be specified.

Currently, the best way to run GELATO is using the wrapper scripts in the in the Convenience subdirectory. The scripts are executable and can be called directly. Ensure that you are in the GELATO conda environment before running any of the scripts. These scripts can be copied to your working directory.

There are two wrappers for GELATO, one for running on a single spectrum, and one for running on a list of spectra.

* "run_GELATO_single.py"

   This script is designed to run GELATO over a single object. This takes 3 positional arguments, the path to the parameters file, the path to the spectrum, and the redshift of the object. You can copy this file into your working directory once GELATO has been installed into the conda environment.

  ```bash
  python run_GELATO_single.py PARAMS.json spectrum.fits 0.5
  ```

* "run_GELATO_multi.py"

   This script is designed to run GELATO over a list of objects. This takes 2 positional arguments, the path to the parameters file, and the path to the list of objects. You can copy this file into your working directory once GELATO has been installed into the conda environment. While it is called multi, it can be run on a file containing a single object.

  ```bash
  python run_GELATO_multi.py PARAMS.json spectra_with_redshifts.csv
  ```

During a GELATO run, rest equivalent widths and plots can be generated depending on what is specified in the parameter file. However, if you opt out of creating them during the run, you can always create them after using The following scripts. These scripts can be copied to the working directory after the installation. Similarly, these scripts are executable and can be called directly.

* "Plotting.py":

  ```bash
  # For a single plot
  python Plot_from_results.py PARAMS.json --Spectrum spectrum.fits --Redshift 0.5
  ```

  ```bash
  # For multiple plots
  python Plot_from_results.py PARAMS.json --ObjectList spectra_with_redshifts.csv
  ```

* "EquivalentWidth.py"

  ```bash
  python EW_from_results.py PARAMS.json --Spectrum spectrum.fits --Redshift 0.5
  ```

  ```bash
  python EW_from_results.py PARAMS.json --Spectrum spectrum.fits --Redshift 0.5
  ```

The concatenated results, with median parameters and standard deviations, for GELATO can also be created directly the results files in the following manner:

* "Concat_from_results.py"

  ```bash
  python Concat_from_results.py PARAMS.json spectra_with_redshifts.csv
  ```

Running the Example
-------------

We provided an example for running GELATO on a few SDSS spectra.This tutorial assumes you start in the Example directory. First we need to activate our GELATO environment.

```bash
conda activate gelato
```

We can then run the code over the whole data set.

```bash
python ../Convenience/run_GELATO_multi.py ExPARAMS.json ExObjList.csv
```

The output from running the example will be put into 'Results/' and can be compared to the results in the 'Comparison/' directory.

While GELATO is designed to be run in this fashion, an IPython notebook is provided in the Example directory. This can also help with how to access GELATO output.

Result File
-------------

The results are presented in a file ending with "-results.fits". It is a multi extension FITS file. Each extension is named based on its contents and can be retrieved in the following manner:
```python
from astrop.io import fits
print(fits.open('example-results.fits').info())
```
The extensions are:

* SUMMARY: It is a binary FITS table containing the summary of the models. It contains the original spectrum without the bad data points (ivar = 0) along with the total model, ssp continuum, power-law continuum, emission-line model generated from the median parameters. The latter two columns will not appear if they are not included in the final fit of the spectrum.
* PARAMS: It is a binary FITS table where each column represents a parameter with a row for each bootstrap. Redshift and dispersion measurements are given in km/s. Flux measurements are dependent on the input units of the spectrum. In addition, the rest amplitude (RAmp) of the Gaussian is also returned as this can be a more reliable way for computing line detection. If calculated, rest equivalent widths are given in Angstroms. If results are concatenated, errors are the standard devations on the recovered parameters from the bootstraps. Coefficients on the SSP continuum models are in the units of the SSP models.


Parameter File
-------------

The behavior of GELATO is controlled entirely by the JSON parameters file. And example parameter file is included in the Example directory.

* Outfolder: This parameter is the path to the output directory.
* VacuumWav: Are the spectra being fit in air or vacuum wavelengths.  
* RandomSeed: The seed used as input to NumPy for random number generation.
* ContinuumRegion: The border around emission lines in velocity space that will be excluded when fitting the continuum initially.
* LineRegion: The border around an emission line in velocity that must be contained within the spectrum in order to be fit (km/s). This region is also used to estimate the initial height of the line.
* NBoot: Number of bootstrap iterations to constrain error on parameters.
* FThresh: F-test threshold to incorporate additional model parameters.
* NProcess: Number of processes to open with python multiprocessing. Set equal to 1 to use only a single thread.
* Plotting: Produce plots or not.
* FlamUnits: String containing the units of the spectrum flux for purposes of plotting, can accept LaTeX syntax.
* CalcEW: To calculate (rest) equivalent widths or not.
* Concatenate: To concatenate the results of a multiple GELATO run or not.
* Overwrite: Overwrite the results of a previous GELATO run.
* Verbose: To print GELATO output.
* EmissionGroups: Dictionary of emission lines to be fit by GELATO. The structure of this dictionary is crucial to the operation of GELATO. The following section details the format of this dictionary.

Emission Line Dictionary
-------------

Here we describe the format of the emission line dictionary.

1. The emission groups dictionary is made up of Groups. All spectral features in the same group can be set to share a common redshift, a common dispersion, or neither.  This means, during fitting, their redshifts or dispersions can be forcibly tied to be equal.

      * Each group has a Name, which controls how its parameters appear in the output.
      * Each group has TieRedshift flag, which controls if the redshifts of all the group elements are tied or not.
      * Each group has TieDispersion flag, which controls if the dispersions of all the group elements are tied or not.
      * Finally, each group is made out of a list of species.

2. Each Group contains a list of Species. All spectral features in the same Species will share a redshift velocity and dispersion. This means, during fitting, their velocity dispersions and redshifts will be forcibly tied to be equal.

      * Each species has a name, which controls how its parameters appear in the output.
      * Each species has a Flag. The integer flag controls will additional parameters GELATO will attempt to add to the spectral features of this species. The value in each bit, from right to left (increasing order of magnitude), is a boolean flag for each kind of additional component, which can be found in the Additional Components section of the README and the "AdditionalComponents.py" file.
      * Each species has a FlagGroup. Each additional component must be associated with a Group, which may or may not be the same Group as the original species. An additional list is passed to each species, specifying where each additional component that will attempt to be added must go. The group must exist, even if empty, as it needs to be created with the flags. (Note: This means the sum of the bits of the flag must be equal to the length of the FlagGroups list.)

3. Each Species contains a list of lines. Lines can be set to have relative fluxes.  This means, during fitting, their fluxes will be tied to have the given relative values.

      * Lines have a Wavelength. This is the rest wavelength of the line (same units as spectrum wavelength).
      * Lines have a RelStrength. This is a relative strength to the other members of the species. If set to null, it will have an independent flux.

The "ExampleParameters.json" file in the Example directory gives a good example of how to take advantage of these features. This is not mean to represent a physically accurate sample, but to give an idea of all of the features of the code. It consists of four groups:

1. Name: AGN. Here we might put our AGN lines which in this case we want to share redshifts but we don't want to share dispersions. It has a list of species:
   1. Name: [SII]. These features will share the same velocity dispersion and redshift. There are no flags on this component, so the list is empty. It is made out of two lines.
      * A line with a rest wavelength of 6716.44 and its flux is left free.
      * A line with a rest wavelength of 6730.82 and its flux is left free. This means the line fluxes are completely independent.
   2. Name: [NII]. These features will share the same velocity dispersion and redshift. These have been flagged with a 2, or in binary 10. This corresponds to an "Outflow" component. This additional component will be placed the group named "Outflow". It is made out of two lines.
      * A line with a rest wavelength of 6583.45 and a relative flux of 1.
      * A line with a rest wavelength of 6548.05 and a relative flux of 0.34. This means this line will always have 0.34/1 times the flux of the first line.
   3. Name: [OIII]. These features will share the same velocity dispersion and redshift. These have been flagged with a 2, or in binary 10. This corresponds to an "Outflow" component. This additional component will be placed the group named "Outflow". It is made out of three lines.
      * A line with a rest wavelength of 5006.84 and a relative flux of 1.
      * A line with a rest wavelength of 4958.91 and a relative flux of 0.35. This means this line will always have 0.35/1 times the flux of the first line.
      * A line with a rest wavelength of 4363.21 and its flux is left free.
2. Name: SF. Here these features are set to share redshifts and dispersions. It has a list of species:
   1. Name: [OI].  There are no flags on this component, so the list is empty. It is made out of two lines.
      * A line with a rest wavelength of 6300.3 and a relative flux of 3.
      * A line with a rest wavelength of 6463.78 and a relative flux of 1. This means this line will always have 1/3 times the flux of the first line.
   2. Name: [OII]. A singlet line with no flags, so the list is empty.
      * A line with a rest wavelength of 3727.43 and its flux is left free.
3. Name: Balmer. Here these features will not share redshifts and dispersions. It has a list of species:
   1. HI. These features will share the same velocity dispersion and redshift. These have been flagged with a 1, which corresponds to a "Broad" component. This additional component will be placed the group named "Balmer". It is made out of one line.
      * A line with a rest wavelength of 6562.79 and its flux is left free.
      * A line with a rest wavelength of 4861.28 and its flux is left free.
      * A line with a rest wavelength of 4340.47 and its flux is left free.
4. Name: Outflow
  If more than one component lands in this group, they will share redshifts and dispersions. E.g. if "Outflow" lines are accepted (from [NII] and [OIII]), they will share the same redshift and dispersion by design.

Here is table showing the hierarchy of the Emission Groups Parameter for the "PARAMS.json" file. A script, params_to_TeX.py, is provided in the Convenience directory can turn a Emission Groups dictionary in a Parameter file into a LaTeX table.

  ```bash
  python params_to_TeX.py PARAMS.json
  ```

![Image of PARAMS](./Images/PARAMS.png)

Models
-------------

* Emission Line Model: Emission lines are modeled as Gaussians parametrized with a redshift, a flux, and a dispersion (in km/s). They are forced to have a positive flux. The default value of the velocity dispersion of the line is set to 150 km/s, while it is bounded between 60 km/s and 1000 km/s. This default can be adjusted in the "CustomModels.py" file.

* Continuum SSP Model: The continuum is modeled as the sum of E-MILES SSP models. In total, 15 SSP models are used to build a continuum. The normalization coefficients are named for each SSP model.

* Continuum Power Law Model: An additional power law continuum is attempted to be fit in addition to the SSP models. It is parametrized with a power law index, a normalization coefficient, and a scale (y = coeff*(x/scale)**(-index)). The power law index has a default value of 1.5. The scale is set as the 20th percentile of the wavelength values for which there is a good flux value. i.e. np.nanpercentile(10**loglam[ivar > 0],20).

Additional Components
-------------

The current supported additional components are:

1. Broad Component: The broad components are modeled as Gaussians. They are forced to have a positive flux. The default value of the velocity dispersion of the line is set to 1200 km/s, while it is bounded between 1000 km/s and 10000 km/s.
2. Outflow Component: They are forced to have a positive flux. The outflow components are modeled as Gaussians. The default value of the velocity dispersion of the line is set to 500 km/s, while it is bounded between 500 km/s and 1000 km/s.

In order to have GELATO attempt to fit an emission line with an additional component, the line must be flagged in the parameters file, described in the section below. The flag is an integer, whose bitwise digits describe if a specific additional component should be tried. Examples for all possible combinations are given in the figure following the description of the EmissionGroups parameter.

GELATO submodules
-------------

* ConstructParams.py

  Routines for turning the PARAMS.json file into a python dictionary, and verifying that it is in the correct format.

* GELATO.py

  Main GELATO function that calls and coordinates the whole operation.

* SpectrumClass.py

  A class the defines how a spectrum is loaded. GELATO was designed for SDSS spectra. To load in any other kind of spectrum, you can edit the way this class is initialized.

* BuildModel.py

  This file handles the construction of models and for tying the various parameters together as outlined by the emission line dictionary.

* CustomModels.py

  Here are where the custom models used in GELATO are defined. Here exists a gaussian emission line model, the SSP continuum, and the power law continuum. The parameters for each model are defined with respect to the rest frame, but the output of the model is in the observed frame. This is where the velocity width limits on emission features can be set.

* AdditionalComponents.py

  Here are where the additional components are defined along with their bit flag positions. In order to add extra additional components, this file can be easily extended to include more models. This is where the velocity dispersion limits on additional components can be modified.

* FittingModel.py

  Here are the scripts for fitting GELATO generated models and for testing the inclusion of additional parameters. To change the fitting algorithm, this file can be edited.

* ModelComparison.py

  Here are scripts for model comparison and selection, including F-tests and AIC calculation.

* Plotting.py

  Here are the scripts for creating and saving figures of the fits. Can also be run directly on GELATO results in order to create figures after the fact. The plots can either be plotted as the components or the final result only.

* EquivalentWidth.py

  Here are the scripts for creating and saving emission line EW. Can also be run directly on GELATO results in order to generate EW after the fact. Rest equivalent widths are generated by assuming a flat continuum at the height of the continuum at the emission line center.

* Concatenate.py

  Scripts for concatenating results from a multi GELATO run. Can also be run independently on results after the fact.

License
-------------

GELATO is an open-source software available under the GNU General Public License 3. In a nutshell, this code can be used and distributed by anyone, but any code that includes GELATO must also be distributed freely and openly (see LICENSE file for details).

FAQ
-------------

**How can I load spectra from other sources?**

*By editing the SpectrumClass.py file, you can customize how spectra are loaded into GELATO. However it might be easier to convert your spectrum to follow the SDSS convention.*

**What are the units?**

*ContinuumRegion and LineRegion are quoted velocity space and are given in km/s. Otherwise, the code is agnostic to the flux units. The wavelengths must be given in Angstroms.*

**Do you mean velocity offsets, not redshifts?**

*Each emission line is characterized by a redshift, which is trivial to convert to a velocity offset once a reference line is chosen. However this requires the user to choose a reference line. GELATO remains agnostic to this procedure and simply returns the redshift of each line.*

**What does it mean if GELATO failed on an object?**

*If you receive the notice that GELATO failed on an object, this means the SVD algorithm at the heart of the SciPy fitting routine failed. In my experience, this issue is operating system, CPU, and linear algebra library dependent. I recommend running the object through GELATO again with a different random seed or on a different machine.*
