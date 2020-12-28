# Netcdf Editor App

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/CEREGE-CL/netcdf_editor_app/main?filepath=app.ipynb) [![Heroku App Status](https://heroku-shields.herokuapp.com/netcdf-editor-app)](https://netcdf-editor-app.herokuapp.com) [![Generic badge](https://img.shields.io/badge/Docs-up-brightgreen.svg)](https://cerege-cl.github.io/netcdf_editor_app/)

Example can be found here: https://netcdf-editor-app.herokuapp.com/

## Introduction

This is a small web based app that allows users to interactively modify NetCDF files. 

Currently 3 methods for modifying values has been implemented:
- __Absolute__: All selected values are set to the replacement value
- __Relatif__: Add replacement value to all selected values
- __Percentage__: Add percentage of each selected value

## App Workflow

1. Load NetCDF File
1. Choose Variable
1. Choose Zones
    - Choose zone on map
    - Choose part of distribution
1. Click Apply
1. Save / Download

## Installation

### Docker

1. Run `docker-compose -f "docker-compose.yml" up  --build`

## Developpement

### testing

The test suite uses pytest and coverage.

To run the tests either run
```shell
python -m pytest tests/
```
or:
```shell
coverage run -m pytest && coverage report
```

