FROM python:3.8

# Setup the enviroment
COPY requirements.*txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt

COPY services/python_worker/receive.py services/python_worker/steps.py /usr/src/app/

COPY netcdf_editor_app /usr/src/app/netcdf_editor_app
WORKDIR /usr/src/app

ENTRYPOINT [ "python", "/usr/src/app/receive.py"]