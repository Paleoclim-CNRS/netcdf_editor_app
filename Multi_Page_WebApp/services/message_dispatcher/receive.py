#!/usr/bin/env python
import pika
import sys
import os
import json

from constants import tasks, invalidates, invalidated, no_params


def send_task(task, body, ch):
    for worker, worker_tasks in tasks.items():
        if task in worker_tasks:
            routing_key = "preprocessing" + "." + task + "." + worker
            ch.basic_publish(
                exchange="preprocessing",
                routing_key=routing_key,
                body=body,
                properties=pika.BasicProperties(
                    delivery_mode=2,  # make message persistent
                ),
            )
            print(
                " [x] Sent message to {} {}".format(routing_key, body),
                flush=True,
            )
            return


def main():
    connection = pika.BlockingConnection(
        pika.ConnectionParameters(host=os.environ["BROKER_HOSTNAME"])
    )

    channel = connection.channel()

    channel.exchange_declare(exchange="preprocessing", exchange_type="topic")
    channel.queue_declare(queue="preprocessing_task_queue", durable=True)

    channel.queue_bind(
        exchange="preprocessing",
        queue="preprocessing_task_queue",
        routing_key="preprocessing.*",
    )
    channel.queue_bind(
        exchange="preprocessing",
        queue="preprocessing_task_queue",
        routing_key="preprocessing.*.done",
    )

    def callback(ch, method, properties, body):
        # print(" [x] ch: ", ch, flush=True)
        # print(" [x] method: ", method, flush=True)
        # print(" [x] properties: ", properties, flush=True)
        print(" [x] Received %r" % body.decode(), flush=True)

        routing_key = method.routing_key
        print(" [x] Received %r" % routing_key, flush=True)

        # task is second object
        task = routing_key.split(".")[1]

        # Deal with new tasks
        if len(routing_key.split(".")) == 2:
            send_task(task, body, ch=ch)

        if len(routing_key.split(".")) == 3:
            if task in invalidates.keys():
                tasks_to_process = None
                if json.loads(body).get("next_step_only", False):
                    tasks_to_process = invalidates[task]
                else:
                    tasks_to_process = invalidated(task)
                print(f"[X] task: {task} invalidates {tasks_to_process}", flush=True)
                _id = json.loads(body)["id"]
                send_task(
                    "invalidate",
                    json.dumps({"id": _id, "tasks": tasks_to_process}),
                    ch=ch,
                )

                for step_invalidated in tasks_to_process:
                    _body = {"id": _id, "invalidated": "yes"}
                    if step_invalidated in no_params:
                        _body["has_params"] = "no"
                    send_task(step_invalidated, body=json.dumps(_body), ch=ch)

        print(" [x] Done", flush=True)
        ch.basic_ack(delivery_tag=method.delivery_tag)

    channel.basic_qos(prefetch_count=1)
    channel.basic_consume(
        queue="preprocessing_task_queue", on_message_callback=callback
    )

    print(" [*] Waiting for messages. To exit press CTRL+C", flush=True)
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
