FROM python:3.8

# Setup the enviroment
RUN pip install pika numpy

COPY services/message_dispatcher/receive.py /tmp/
COPY netcdf_editor_app/constants.py /tmp/

ENTRYPOINT [ "python", "/tmp/receive.py"]