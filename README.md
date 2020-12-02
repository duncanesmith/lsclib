# lsclib

## Summary 
Luminescent solar concentrators (LSCs) enhance the power output of solar cells via luminescent emission and internal reflection.
LSCs have long been speculated as BIPV due to their innate architectural flexibility and vast potential for improvement in PV efficiency.
However, it is more difficult to model LSCs as compared with solar panels, and this has limited their integration commercially. lsclib
is a python-based repository hosted on GitHub that employs the Monte Carlo ray-tracing method of radiative transport to effectively model LSCs.

lsclib hopes to short-circuit the learning curve associated with breaking into the field, and present results in both an academic and
commercial context. This repository will continue to become more sophisticated, but for now relies heavily upon the paper submitted for publishing
entitled "An open-source Monte Carlo ray-trace simulation tool for luminescent solar concentrators with validation studies employing scattering phosphor films".
Install the [lsclib package](https://pypi.org/project/lsclib/) to get started!

While lsclib is growing and improving, the visitor to this repository is also encouraged to visit [pvtrace](https://github.com/danieljfarrell/pvtrace). 
lsclib has many advantages, but pvtrace is quite extensive and may be a better fit for certain applications.
	
## Requirements
To ensure that lsclib can run properly, the first thing you'll want to do is ensure you have installed the required packages seen below.
- [matplotlib](https://pypi.org/project/matplotlib/)
- [numpy](https://pypi.org/project/numpy/)
- [pandas](https://pypi.org/project/pandas/)
- [scipy](https://pypi.org/project/scipy/)
- [Shapely](https://pypi.org/project/Shapely/)

## Getting Started
Before you start running any code, you should take a look at the following python files: "run", "lsc_classes", and "lsc_calcs". The "run" module currently houses
one example function called "wedge" which can produce the data seen in the referenced paper. "lsc_classes" holds the python Classes that make up an LSC, and you'll
see that the "lsc_classes" module is referenced frequently by the "run" module. "lsc_calcs" houses the majority of theoretical calculations that the model relys upon.
You'll see that the "lsc_calcs" module is referenced frequently by the "lsc_classes" module. These and all other modules used are commented - open the lsclib folder
to view everything.

A good line of code to get you started (once everything is installed and imported) is: "LSC = lsclib.run.wedge(1000)", which will create an LSC object with the affiliated
attributes produced for 1000 trials under normal insolation. For example, "LSC.Isc_cell" will report the resulting short-circuit current of the solar cell within an LSC.