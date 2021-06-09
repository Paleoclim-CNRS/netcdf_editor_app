import pika
from flask import g
import os
import json

from datetime import datetime


def get_connection():
    if "connection" not in g:
        g.connection = pika.BlockingConnection(
            pika.ConnectionParameters(os.environ["BROKER_HOSTNAME"])
        )

    return g.connection


def close_connection(e=None):
    connection = g.pop("connection", None)

    if connection is not None:
        connection.close()


def init_app(app):
    app.teardown_appcontext(close_connection)


def get_channel():
    connection = get_connection()
    return connection.channel()


def send_preprocessing_message(routing_key, message=None):

    channel = get_channel()

    channel.exchange_declare(exchange="preprocessing", exchange_type="topic")

    channel.basic_publish(
        exchange="preprocessing",
        routing_key="preprocessing." + routing_key,
        body=json.dumps(message),
        properties=pika.BasicProperties(
            delivery_mode=2,  # make message persistent
        ),
    )
    print(f" [x] {datetime.now()} Message Sent {json.dumps(message)}", flush=True)