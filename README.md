# lsclib

[![DOI](https://zenodo.org/badge/313684595.svg)](https://zenodo.org/badge/latestdoi/313684595)

## Summary 
Luminescent solar concentrators (LSCs) enhance the power output of solar cells via luminescent emission and internal reflection.
LSCs have long been speculated as BIPV due to their innate architectural flexibility and vast potential for improvement in PV efficiency.
However, it is more difficult to model LSCs as compared with solar panels, and this has limited their integration commercially. lsclib
is a python-based repository hosted on GitHub that employs the Monte Carlo ray-tracing method of radiative transport to effectively model LSCs.
Additionally, the Monte Carlo code found here can be used to generate inputs for LSC energy estimates.

lsclib hopes to short-circuit the learning curve associated with breaking into the field, and present results in both an academic and
commercial context. This repository will continue to become more sophisticated, but for now relies heavily upon the paper entitled: 
"An Open-source Monte Carlo Ray-Tracing Simulation Tool for Luminescent Solar Concentrators With Validation Studies Employing Scattering Phosphor Films".
This paper can be [found on MDPI](https://www.mdpi.com/1996-1073/14/2/455).
Download the repository to get started!

While lsclib is growing and improving, the visitor to this repository is also encouraged to visit [pvtrace](https://github.com/danieljfarrell/pvtrace). 
lsclib has many advantages, but pvtrace is quite extensive and may be a better fit for certain applications. One primary distinction, however, is that the Monte Carlo tool found here will generate data necessary to compute annual energy estimates for LSCs, relying on functions found primarily in [pvlib python](https://pvlib-python.readthedocs.io/en/stable/). The recently published paper in [Renewable Energy](https://www.sciencedirect.com/science/article/abs/pii/S0960148121018759?via%3Dihub), "Towards a Standard Approach for Annual Energy Production of Concentrator-based Building-integrated Photovoltaics" describes the approach used to run energy estimates relying upon pvlib.
	
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

Once you've downloaded the files from GitHub, navigate to the "run.py" file and run it. Then you can type "LSC = wedge(1000)", which will create an LSC object with the affiliated
attributes produced for 1000 trials under normal insolation. For example, "LSC.Isc_cell" will report the resulting short-circuit current of the solar cell within an LSC. The "energy_estimate.py" file will run energy estimates for both LSCs (forecast_lsc) and solar panels (forecast_solar_panel) as well. Existing csv files for LSC optical efficiency and spectral mismatch factor can be leveraged, or recreated by running a matrix of results through the Monte Carlo code.
