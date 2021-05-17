FROM python:3.8

# Setup the enviroment
COPY requirements.*txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt

COPY services/python_worker/receive.py services/python_worker/steps.py /usr/src/app/

COPY flask_app /usr/src/app/flask_app
RUN pip install /usr/src/app/flask_app
COPY python_tools /usr/src/app/python_tools
RUN pip install /usr/src/app/python_tools
WORKDIR /usr/src/app

ENTRYPOINT [ "python", "/usr/src/app/receive.py"]