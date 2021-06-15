#!/usr/bin/env python
from datetime import datetime
import pika
import os
import sys
import steps  # noqa: F401
import json

from climate_simulation_platform.db import step_parameters, save_step, step_seen
from climate_simulation_platform import create_app


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


def main():
    connection = pika.BlockingConnection(
        pika.ConnectionParameters(host=os.environ["BROKER_HOSTNAME"])
    )

    app = create_app()

    channel = connection.channel()

    channel.exchange_declare(exchange="preprocessing", exchange_type="topic")
    channel.queue_declare(queue="preprocessing_python_task_queue", durable=True)

    channel.queue_bind(
        exchange="preprocessing",
        queue="preprocessing_python_task_queue",
        routing_key="preprocessing.*.python",
    )

    def callback(ch, method, properties, body):

        routing_key = method.routing_key
        print(
            f" [x] {datetime.now()} Received message from {routing_key} with body: {body.decode()}",
            flush=True,
        )
        func = routing_key.split(".")[1]
        body = json.loads(body.decode())
        params = func_params(func, body)
        print(f"{datetime.now()} Params: {params}", flush=True)
        if params is not None:
            _id = body["id"]
            if func != "invalidate":
                with app.app_context():
                    save_step(_id, func, params, up_to_date=False)
            eval(f"steps.{func}({params})")
            if func != "invalidate":
                with app.app_context():
                    save_step(_id, func, params, up_to_date=True)

            routing_key_done = ".".join([*routing_key.split(".")[:2], "done"])
            channel.basic_publish(
                exchange="preprocessing",
                routing_key=routing_key_done,
                body=json.dumps(body),
                properties=pika.BasicProperties(
                    delivery_mode=2,  # make message persistent
                ),
            )
            print(
                " [x] {} Sent message to {} {}".format(
                    datetime.now(), routing_key_done, body
                ),
                flush=True,
            )
        print(f" [x] {datetime.now()} Done", flush=True)
        ch.basic_ack(delivery_tag=method.delivery_tag)

    channel.basic_qos(prefetch_count=1)
    channel.basic_consume(
        queue="preprocessing_python_task_queue", on_message_callback=callback
    )

    print(
        f" [*] {datetime.now()} Waiting for messages. To exit press CTRL+C", flush=True
    )
    channel.start_consuming()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("Interrupted")
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
