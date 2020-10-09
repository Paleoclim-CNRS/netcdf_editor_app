FROM jupyter/minimal-notebook

COPY requirements.txt /tmp/
RUN pip install -r /tmp/requirements.txt

WORKDIR /app
COPY app.ipynb .

CMD ["panel", "serve", "app.ipynb"]
