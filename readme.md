# solarpy - Solar Cell and Radiation Degradation Model and Prediction

The solarpy package includes both EQFLUX and DDD degradation models.  Using the functions in this library an user can either perform such functions as fitting radiation degradation data, generate relative damage coefficients (RDCs), derive DDD fitting coeffecients as well as other functions.  Included in this package is the EQFLUX fully translated from the original Fortran code in The GaAs Radiation Handbook. Also included is the radiation degradation data of the ASEC GaAs solar cells. The data was extracted from teh plots in the handbook and is available for testing and examples.

## Getting Started

1.) Once you have download the .zip file and extracted the files, open up the Anaconda prompt. <br />
<br />
2.) Install using the command "pip install 'path to package'". Using pip to install will install all dependencies <br />
2a.) Alternatively if you plan to actively develop pearl, you can install using "pip install -e" or cd into the solarpy directory to the setup.py and use "python setup.py develop" (for advanced users). <br />
<br />
3.) The solarpy library will automatically load all dependencies if you use pip.<br />
<br />
4.) Check to see if solarpy is installed by typing "python" in the command line. Next type "import solarpy as sy" and type "sy.GaAs_Electron_RDC". You should see the print out for the GaAs Electron Rdc from The GaAs Radiation Handbook.<br />
<br />
6.) If solarpy is correctly installed you are good to go...if not email don.walker@aero.org

### Prerequisites

python 3.10

### solarpy DOCUMENTATION

Docs can be found here: <a href = "https://the-aerospace-corporation.github.io/solarpy/" target="_blank">Solarpy Docs.</a>
```
numpy
scipy
pandasgit
matplotlib
julian
profilehooks
```

