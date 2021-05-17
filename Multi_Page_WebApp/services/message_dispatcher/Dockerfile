FROM python:3.8

# Setup the enviroment
RUN pip install pika numpy

COPY services/message_dispatcher/receive.py /tmp/
COPY flask_app/climate_simulation_platform/constants.py /tmp/

ENTRYPOINT [ "python", "/tmp/receive.py"]