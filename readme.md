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

python 3.8

### solarpy DOCUMENTATION
https://aerosource2.aero.org/bitbucket/pages/SOLARPY/solarpy/master/browse/docs/_build/html/index.html

```
numpy
scipy
pandas
matplotlib
julian
profilehooks
```

### Making Changes to the Project
Use the following basic git commands to clone, checkout your new branch, and push your changes.

#### Cloning
```
git clone https://DW28596@aerosource2.aero.org/bitbucket/scm/solarpy/solarpy.git
```

#### Switching Branches
```
git checkout [branch name]
```

#### Adding Changes to git
```
// Add the files you with to submit changes for
git add [file name/path]

// Save Changes on Local Machine
git commit -m "My Change Log Message"

// Sync Changes with Remote Machine
git push
```

#### Merge changes with *main* branch
1. Open a pull request to merge your branch with dev branch
2. Add reviewers
3. Add any additional message for pull request


### Installing
Use pip to install solarpy.  Using pip will download all required packages.

Example
pip install 'c:\Users\solarpy'

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags).

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
