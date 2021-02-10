FROM python:3.8

# Setup the enviroment
COPY requirements.*txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt

COPY netcdf_editor_app /usr/src/app/netcdf_editor_app
WORKDIR /usr/src/app

ENV BOKEH_RESOURCES=cdn 
ENV BOKEH_ALLOW_WS_ORIGIN=localhost:9090

ENTRYPOINT [  "python", "-m", "panel",  "serve",  "--port=5006", "--allow-websocket-origin='*'", "--websocket-max-message-size=100000000", "netcdf_editor_app/holoviews/internal_oceans.py", "netcdf_editor_app/holoviews/value_changer.py", "netcdf_editor_app/holoviews/passage_problems.py" ]