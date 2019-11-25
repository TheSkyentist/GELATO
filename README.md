HQUAILS
========
*Handy QUAsar emissIon Line fitS (pronounced Quails) by Raphael Hviding*
-------------

HQUAILS is a Python code designed to fit emission lines in the spectra of active galactic nuclei. In particular, it was built in order to fit AGN spectra where many of the parameters of the emission lines are tied with respect to one another. HQUAILS attempts to automate this process. 

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

License
-------------
HQUAILS is an open-source software available under the GNU General Public Licence. In a nutshell, this code can be used and distributed by anyone, but any code that includes HQUAILS must also be distributed freely and openly (see LICENCE file for details).