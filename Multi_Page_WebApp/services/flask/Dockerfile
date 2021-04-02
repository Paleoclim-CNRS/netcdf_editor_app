FROM python:3.8

# Setup the enviroment
COPY requirements.*txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt
RUN pip install gunicorn

COPY netcdf_editor_app /usr/src/app/netcdf_editor_app
ENV FLASK_APP=/usr/src/app/netcdf_editor_app

WORKDIR /usr/src/app
RUN python -m flask init-db

# ENTRYPOINT [ "waitress-serve", "--call", "--listen=${FLASK_RUN_HOST}:${FLASK_RUN_PORT}", "'netcdf_editor_app:create_app'" ]
# ENTRYPOINT [ "waitress-serve", "--listen=*:5000",  "--call", "netcdf_editor_app:create_app" ]