## Introduction
{:.no_toc}
The NetCDF Editor App is maintained by [CEREGE-CL](https://github.com/CEREGE-CL). The goal of the app is to be able to interactively visualize and adjust netcdf files for use for deep time simulations.

## Table of Contents
{:.no_toc}
1. TOC
{:toc}

## Deployements

### Local
The repository is setup to be able to run on local hardware (user's laptop or on premise cloud). To run locally simply clone this repository and run:

```docker-compose -f "docker-compose.yml" up --build```

### Cloud

# IMAGE GOES HERE

#### [Heroku](https://netcdf-editor-app.herokuapp.com)

[![Heroku App Status](https://heroku-shields.herokuapp.com/netcdf-editor-app)](https://netcdf-editor-app.herokuapp.com)

Heroku is a free service that is automatically deployed from the github repo when the `main` branch is updated. 

Being a free service ressource are limited and the app is laggy. This is maybe the best way to get a gist of what is happening but not useful for carrying out work. 

You can test the app [here](https://netcdf-editor-app.herokuapp.com)

#### [Google Cloud](https://netcdf-editor.ew.r.appspot.com/app)

We are testing using Google Cloud (Google App Engine) to build and run the App. In the same manner a trigger has been setup so that when `main` is updated a build and deployment occurs automatically.

Currently the App is running on free credits, however when these credits run out we have to decide if we carry on using this service or not.

You can test the app [here](https://netcdf-editor.ew.r.appspot.com/app)

## Interface

### Browse

The entrypoint to the app is to upload a file using the available upload tool, seen below.

![](img/browse.png)

> Currently Google Cloud only accept files less than roughly 20mb, however the limit for local deployement has been set to 100mb and is modifiable.

Under the hood the import and manipulation of the netCDF files is done by [xarray](https://xarray.pydata.org/en/stable/). 

```python
import xarray as xr
data = xr.open_dataset(FILENAME)
```

#### Curvilinear Coordinates

Curvilinear coordinates are supported __HOWEVER__ they are reprojected into _ij_ space. This is done because calculating selections on a non uniform grid is costly. Below the curvilinear grid is on the left, notice the areas in black that were not initially covered.

![](img/curvilinear_import.png)
