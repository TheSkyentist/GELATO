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

In your working directory, you need to copy the "PARAMS.py" and "matplotlibrc" file. The former controls the HQUAILS settings, while the latter controls the plotting settings.

Running HQUAILS
-------------

How it works
-------------



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