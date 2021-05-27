#!/usr/bin/env python
import pika
import os
import steps  # noqa: F401
import json

from climate_simulation_platform.db import step_parameters, save_step, step_seen
from climate_simulation_platform import create_app


# -*- coding: utf-8 -*-
# pylint: disable=C0111,C0103,R0205

import functools
import logging
import threading

LOG_FORMAT = (
    "%(levelname) -10s %(asctime)s %(name) -30s %(funcName) "
    "-35s %(lineno) -5d: %(message)s"
)
LOGGER = logging.getLogger(__name__)

logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)


def func_params(func, body):
    # If invalidated isn't in keys then this is a "root" call meaning it should be run
    if "invalidated" not in body.keys():
        return body
    # If 'invalidated': 'y(es)' in the body then this means the step has been invalidated
    # It should be rerun IF it has already been run before OR has no params
    # We will rerun it with the same parameters
    if "invalidated" in body.keys() and body["invalidated"].lower() in ["yes", "y"]:
        if "has_params" in body.keys() and body["has_params"].lower() in ["no", "n"]:
            return body
        app = create_app()
        with app.app_context():
            if step_seen(body["id"], func):
                return step_parameters(body["id"], func)
    return None


def run_function(method, body):
    app = create_app()
    routing_key = method.routing_key
    print(" [x] Received %r" % routing_key, flush=True)
    func = routing_key.split(".")[1]
    body = json.loads(body.decode())
    params = func_params(func, body)
    if params is not None:
        _id = body["id"]
        if func != "invalidate":
            with app.app_context():
                save_step(_id, func, params, up_to_date=False)
        eval(f"steps.{func}({params})")
        if func != "invalidate":
            with app.app_context():
                save_step(_id, func, params, up_to_date=True)


def send_response_done(method, body, channel=None):
    if channel is None:
        channel = setup_channel()

    routing_key = method.routing_key
    routing_key_done = ".".join([*routing_key.split(".")[:2], "done"])
    channel.basic_publish(
        exchange="preprocessing",
        routing_key=routing_key_done,
        body=body,
        properties=pika.BasicProperties(
            delivery_mode=2,  # make message persistent
        ),
    )
    print(
        " [x] Sent message to {} {}".format(routing_key_done, body),
        flush=True,
    )


def setup_channel():
    connection = pika.BlockingConnection(
        pika.ConnectionParameters(host=os.environ["BROKER_HOSTNAME"])
    )
    channel = connection.channel()

    channel.exchange_declare(exchange="preprocessing", exchange_type="topic")
    channel.queue_declare(queue="preprocessing_mosaic_task_queue", durable=True)

    channel.queue_bind(
        exchange="preprocessing",
        queue="preprocessing_mosaic_task_queue",
        routing_key="preprocessing.*.mosaic",
    )
    return channel


# def main():
#     # connection = pika.BlockingConnection(
#     #     pika.ConnectionParameters(host=os.environ["BROKER_HOSTNAME"])
#     # )

#     # channel = connection.channel()

#     # channel.exchange_declare(exchange="preprocessing", exchange_type="topic")
#     # channel.queue_declare(queue="preprocessing_mosaic_task_queue", durable=True)

#     # channel.queue_bind(
#     #     exchange="preprocessing",
#     #     queue="preprocessing_mosaic_task_queue",
#     #     routing_key="preprocessing.*.mosaic",
#     # )

#     def callback(ch, method, properties, body):
#         # Acknowldeg we have received the message -> We
#         # do this first because if the process is too long running
#         # Then the broker kicks us out
#         ch.basic_ack(delivery_tag=method.delivery_tag)
#         # print(" [x] ch: ", ch, flush=True)
#         # print(" [x] method: ", method, flush=True)
#         # print(" [x] properties: ", properties, flush=True)
#         print(" [x] Received %r" % body.decode(), flush=True)


#         print(" [x] Done", flush=True)

#     channel.basic_qos(prefetch_count=1)
#     channel.basic_consume(
#         queue="preprocessing_mosaic_task_queue", on_message_callback=callback
#     )

#     print(" [*] Waiting for messages. To exit press CTRL+C", flush=True)
#     channel.start_consuming()


# if __name__ == "__main__":
#     try:
#         main()
#     except KeyboardInterrupt:
#         print("Interrupted")
#         try:
#             sys.exit(0)
#         except SystemExit:
#             os._exit(0)


def ack_message(ch, delivery_tag, method, body):
    """Note that `ch` must be the same pika channel instance via which
    the message being ACKed was retrieved (AMQP protocol constraint).
    """
    if ch.is_open:
        ch.basic_ack(delivery_tag)
        send_response_done(method, body, channel=ch)
    else:
        # Channel is already closed, so we can't ACK this message;
        # log and/or do something that makes sense for your app in this case.
        send_response_done(method, body)


def do_work(conn, ch, method, delivery_tag, body):
    thread_id = threading.get_ident()
    LOGGER.info(
        "Thread id: %s Delivery tag: %s Message body: %s", thread_id, delivery_tag, body
    )
    run_function(method, body)
    cb = functools.partial(ack_message, ch, delivery_tag, method, body)
    conn.add_callback_threadsafe(cb)


def on_message(ch, method_frame, _header_frame, body, args):
    (conn, thrds) = args
    print(" [x] Received %r" % body.decode(), flush=True)
    delivery_tag = method_frame.delivery_tag
    t = threading.Thread(
        target=do_work, args=(conn, ch, method_frame, delivery_tag, body)
    )
    t.start()
    thrds.append(t)


# # Note: sending a short heartbeat to prove that heartbeats are still
# # sent even though the worker simulates long-running work
parameters = pika.ConnectionParameters(host=os.environ["BROKER_HOSTNAME"], heartbeat=5)
connection = pika.BlockingConnection(parameters)


channel = connection.channel()
channel.exchange_declare(exchange="preprocessing", exchange_type="topic")
channel.queue_declare(queue="preprocessing_mosaic_task_queue", durable=True)

channel.queue_bind(
    exchange="preprocessing",
    queue="preprocessing_mosaic_task_queue",
    routing_key="preprocessing.*.mosaic",
)
# Note: prefetch is set to 1 here as an example only and to keep the number of threads created
# to a reasonable amount. In production you will want to test with different prefetch values
# to find which one provides the best performance and usability for your solution
channel.basic_qos(prefetch_count=1)

threads = []
on_message_callback = functools.partial(on_message, args=(connection, threads))
channel.basic_consume("preprocessing_mosaic_task_queue", on_message_callback)

try:
    channel.start_consuming()
except KeyboardInterrupt:
    channel.stop_consuming()

# Wait for all to complete
for thread in threads:
    thread.join()

connection.close()
