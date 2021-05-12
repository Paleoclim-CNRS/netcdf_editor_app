# Netcdf Editor App
[![Generic badge](https://img.shields.io/badge/Docs-up-green.svg)](https://cerege-cl.github.io/netcdf_editor_app/)
[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/netcdf_editor_app/community)



- Single Page netcdf editor app can be found here: https://netcdf.osupytheas.com/
- Multi Page Climate simulation tool can be found here: https://climate_sim.osupytheas.com/

## Introduction

Checkout the [Documentation](https://cerege-cl.github.io/netcdf_editor_app/) for more info.

## Deploying to local server

Images are automatically built using the github actions pipeline see https://github.com/CEREGE-CL/netcdf_editor_app/blob/main/.github/workflows/docker-image.yml these are automatically pushed to dockerhub and are therefore accesible for everyone: https://hub.docker.com/orgs/ceregecl/repositories

To deploy these locally you need to:
1. Copy contents of docker-compose.yaml
1. Replace `${NGINX_PORT}` by the desired port on the server / local machine
1. In `flask_app` and `panel_app` remove the env_file line and add the contents of `config/*.prod` into each corresponding `enviroment`
1. run `docker-compose up --build -d` (you can scale workers with `--scale python_worker=3` where python_worker is the name of the service in docker-compose)

Multipage stack images:
- `ceregecl/netcdf_editor_python_worker`
- `ceregecl/netcdf_editor_flask_app`
- `ceregecl/netcdf_editor_panel_app`
- `ceregecl/netcdf_editor_message_dispatcher`
- `ceregecl/netcdf_editor_nginx`

Single page app:
- `ceregecl/netcdf_editor`


## testing

The test suite uses pytest and coverage.

To run the tests either run
```shell
python -m pytest tests/
```
or:
```shell
coverage run -m pytest && coverage report
```

