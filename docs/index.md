### Introduction
The **Netcdf Editor App** is maintained by [PaleoClimate](https://github.com/Paleoclim-CNRS) team. 

it is composed of 2 underlying applications :

- **Single Page Web App** - for editing netcdf files
- **Multi Page Web App** - for preparing boundary conditions for CM5A2 Earth System Model used for deep time climate simulations.

You can find more informations on the application architecture [here](https://paleoclim-cnrs.github.io/documentation-processing-app/Netcdf_Editor_App_under_the_hood/).

There have already been a succession of iterations of the tool, starting out as a [Single Page](/netcdf_editor_app/single) and evolving into a more complex [Multi Page](/netcdf_editor_app/multi).

<div class='alert alert-info'>
    If you are unsure which tool to use the Multipage web app does everything the single page web app does and is under developpement. It is slightly more ressource heavy (has more tools) but is just as easy to install.
</div>


### Single Page Web App
To be able to quickly manipulate a file on the fly the single page web app is the preferred method.

The source code can be found [here](https://github.com/Paleoclim-CNRS/netcdf_editor_app/tree/main/Single_Page_WebApp).

### Multi Page Web App
However the **Multi Page** allows manipulating files on the fly as well as a number of other tasks:
- Regridding
- Run Routing Code to generate routing files, Bathymetry file (Paleorca), ...
- Generate PFT files
- Generate AHMCoef and Geothermal Heatflow

The Multipage web is built using flask and stores your files in a local database meaning you can come back to your work later.

Also the Multi page app is built with scalability and flexibility in mind, for example it is _trivial_ to add routines to the interface that use multiple languages (Python, Fortran, ...)

This tool also has aspirations of being the turn table in simulation workflows allowing not only the creation of boundary conditions but also more broader tasks such as post processing or simulation lookups.

The source code can be found [here](https://github.com/Paleoclim-CNRS/netcdf_editor_app).

### Installation
You can find documentation on how to deploy these applications locally [here](https://paleoclim-cnrs.github.io/documentation-processing-app/Netcdf_Editor_App_deploy_locally/).