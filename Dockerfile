FROM jupyter/minimal-notebook

COPY requirements.txt /tmp/
RUN pip install -r /tmp/requirements.txt

WORKDIR /app
COPY app.ipynb .

EXPOSE 8080
ENTRYPOINT ["panel", "serve", "--port", "8080", "app.ipynb"]
