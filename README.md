# Introduction to the **ocean** GitHub repository

This [repo](https://github.com/robfatland/ocean) is a collection of Python Jupyter notebooks and some supporting datasets. 
The objective is to create a learning resource for expanding on sensor data visualizations
and to 'lean in' to a synoptic view that combines or compares with data from other sources.


[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/robfatland/ocean/HEAD)


This synoptic ocean data repository is centered on the Shallow Profilers of the
Regional Cabled Array within the Oceans Observing Initiative. The data extend to 
ARGO floats, remote sensing sea surface data and GLODAP. 


Clicking the binder link above creates a sandbox instance of this Jupyter notebook collection. 
Each notebook indicates whether it "runs" (executes Python code in a meaningful way) in binder.
However even when this is not the case: The static cotent may still be of interest (plots of 
data and so forth). 


The data time range centers on January 2019 in order to keep the datasets small; so that they may
in turn be included in the GitHub repository within data volume limits. Examples of datasets too
large to fit include spectrophotometers, hydrophones, ADCP, GLODAP and optical instruments.


The core data are collected in the photic zone by a "shallow profiler" that cycles through the 
upper 200 meters of the water column nine times per day. These nine excursions include two 
intervallic ('with pauses') profiles that accommodate equilibration of the pH sensor. 


### Purpose

This repository illustrates ocean data science. The code is Python broken up as topic
notebooks. Data courtesy of 
[OOI Regional Cabled Array](https://interactiveoceans.washington.edu), NASA, ARGO and other research programs.
To explore: Click on the binder badge above to launch a sandbox copy.
Typical start time will be about 3 minutes. Once binder
finishes try running the first notebook, ***Ocean 01 A Photic Zone***.
When you are done just close the browser tab and the sandbox will evaporate.




### How was this repository binder-ized?

Three steps

- Reduce the source data for the demo notebook down to a few MB so it "lives" in the repo folder
- Add the subfolder `binder` and the `environment.yml` file given below
- Add the binder badge linking to this repo

Badge code: 

```
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/robfatland/ocean/HEAD)
```


`binder/environment.yml` file: 


```
channels:
  - conda-forge
dependencies:
  - python=3
  - numpy
  - pandas
  - matplotlib
  - netcdf4
  - xarray
  - ffmpeg
```


### Other notes on configuration

* On my PC running an Ubuntu `ssh` shell I was obliged to run `conda install dask`.

