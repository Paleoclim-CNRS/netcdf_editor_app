FROM python:3.8

RUN apt-get update && apt-get install -y libgeos-dev

COPY requirements.*txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt

WORKDIR /usr/src/app
COPY app.ipynb .

ENTRYPOINT python -m panel serve --num-procs=1 --port=8080 --address=0.0.0.0 --allow-websocket-origin=* --websocket-max-message-size=100000000 app.ipynb
