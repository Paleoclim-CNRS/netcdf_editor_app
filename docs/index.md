## Introduction
The NetCDF Editor App is maintained by [PaleoClimate](https://github.com/Paleoclim-CNRS) team. The goal of the app is to be able to interactively visualize and adjust netcdf files for use for deep time simulations, in particular to setup the boundary condtions for climate simulations models.

There have already been a succession of iterations of the tool, starting out as a [single page web app](/netcdf_editor_app/single) and evolving into a more complex [multi page web based tool](/netcdf_editor_app/multi).

<div class='alert alert-info'>
    If you are unsure which tool to use the Multipage web app does everything the single page web app does and is under developpement. It is slightly more ressource heavy (has more tools) but is just as easy to install.
</div>


#### Single Page
To be able to quickly manipulate a one file on the fly the single page web app is the preferred method see the [installation instructions](/netcdf_editor_app/single#deployements).

The source code can be found [here](https://github.com/Paleoclim-CNRS/netcdf_editor_app/tree/main/Single_Page_WebApp).

<div class='alert alert-warning'>
There is no current developpement on this
</div>

#### Multi Page
However the multipage web app allows manipulating files on the fly as well as a number of other tasks:
- Regridding
- Run Routing Code to generate routing files, Bathymetry file (Paleorca), ...
- Generate PFT files
- Generate AHMCoef and Geothermal Heatflow

The Multipage web is built using flask and stores your files in a local database meaning you can come back to your work later.

Also the Multi page app is built with scalability and flexibility in mind, for example it is _trivial_ to add routines to the interface that use multiple languages (Python, Fortran, ...)

This tool also has aspirations of being the turn table in simulation workflows allowing not only the creation of boundary conditions but also more broader tasks such as post processing or simulation lookups.

The source code can be found [here](https://github.com/Paleoclim-CNRS/netcdf_editor_app).

See [installation instructions](/netcdf_editor_app/multi#installation) for using the app.
