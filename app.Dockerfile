FROM python:3.8

# Setup the enviroment
COPY requirements.*txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt

COPY netcdf_editor_app /usr/src/app/netcdf_editor_app
ENV FLASK_APP=/usr/src/app/netcdf_editor_app

WORKDIR /usr/src/app/netcdf_editor_app
RUN python -m flask init-db

ENTRYPOINT [ "flask", "run" ]