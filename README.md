# Introduction to the **ocean** GitHub repository

This [repo](https://github.com/robfatland/ocean) is a collection of Python Jupyter notebooks with supporting data. 
It is a learning resource for expanding on sensor data visualization in oceanography. As such there are two
driving emphases: Good dataset construction, and comparison of data from different types of sensors.
The latter is referred to herein as a 'synoptic view' of the water column.


Click the following link -- it may take a few minutes to initialize -- to create a *'binder'* sandbox version
of this repository. Once the temporary copy becomes available in your browser, double-click
the available **BioOptics** Jupyter notebook. Additional work resides in the **`Notebooks`** sub-folder.
You may experiment as you like with the binder instance: It is temporary.


[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/robfatland/ocean/HEAD)


The work here centers on **shallow profilers** located in the
[Regional Cabled Array (RCA)](https://interactiveoceans.washington.edu); 
in turn a large component of NSF's Ocean Observatories Initiative (OOI). 
Additional datasets are featured that originate in other programs 
([RCA cruise data](https://alfresco.oceanobservatories.org/), 
ARGO, 
MODIS and so on).


An RCA shallow profiler rests at a depth of 200 meters below the surface and transits the upper 
water column nine times per day. These nine excursions include two 
intervallic ('with pauses') profiles that accommodate equilibration of a pH sensor. 




### How to create the binder feature for a GitHub repository


- Reduce the source data to a few Megabytes so it can reside in the repository
- Add a subfolder `binder` and an `environment.yml` file as shown given below
- Add a binder badge linking to this repo

Binder link: 

```
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/robfatland/ocean/HEAD)
```


Within the `binder` folder the `environment.yml` file reflects necessary libraries / packages: 


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


### Configuration notes

* One may be obliged to `conda install dask`

