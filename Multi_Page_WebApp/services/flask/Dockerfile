FROM python:3.8

# Setup the enviroment
COPY requirements.*txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt
RUN pip install gunicorn

COPY flask_app /usr/src/app/flask_app
RUN pip install /usr/src/app/flask_app
ENV FLASK_APP=/usr/src/app/flask_app/climate_simulation_platform

WORKDIR /usr/src/app
RUN python -m flask init-db

# ENTRYPOINT [ "waitress-serve", "--call", "--listen=${FLASK_RUN_HOST}:${FLASK_RUN_PORT}", "'climate_simulation_platform:create_app'" ]
# ENTRYPOINT [ "waitress-serve", "--listen=*:5000",  "--call", "climate_simulation_platform:create_app" ]