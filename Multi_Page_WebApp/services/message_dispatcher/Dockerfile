FROM python:3.8

# Setup the enviroment
RUN pip install pika

COPY receive.py /tmp/

ENTRYPOINT [ "python", "/tmp/receive.py"]