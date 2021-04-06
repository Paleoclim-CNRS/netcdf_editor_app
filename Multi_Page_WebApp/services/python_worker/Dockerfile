FROM python:3.8

# Setup the enviroment
COPY requirements.*txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt

COPY services/python_worker/receive.py /tmp/
ENTRYPOINT [ "python", "/tmp/receive.py"]