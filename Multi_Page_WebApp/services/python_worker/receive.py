#!/usr/bin/env python
import pika
import os
import sys
import steps
import json

def should_process(func, body):
    if 'invalidated' not in body.keys():
        return True
    if 'invalidated' in body.keys() and body['invalidated'].lower() in ['yes', 'y']:
        # Test to see if task has already been run
        pass
    return False


def main():
    connection = pika.BlockingConnection(
        pika.ConnectionParameters(host=os.environ["BROKER_HOSTNAME"])
    )

    channel = connection.channel()

    channel.exchange_declare(exchange="preprocessing", exchange_type="topic")
    channel.queue_declare(queue="preprocessing_python_task_queue", durable=True)

    channel.queue_bind(
        exchange="preprocessing",
        queue="preprocessing_python_task_queue",
        routing_key="preprocessing.*.python",
    )

    def callback(ch, method, properties, body):
        print(" [x] ch: ", ch, flush=True)
        print(" [x] method: ", method, flush=True)
        print(" [x] properties: ", properties, flush=True)
        print(" [x] Received %r" % body.decode(), flush=True)

        routing_key = method.routing_key
        print(" [x] Received %r" % routing_key, flush=True)
        func = routing_key.split(".")[1]
        body = json.loads(body.decode())
        print(body, type(body), flush=True)
        if should_process(func, body):
            eval(f"steps.{func}({body})")

            routing_key_done = ".".join([*routing_key.split(".")[:2], "done"])
            channel.basic_publish(
                exchange="preprocessing",
                routing_key=routing_key_done,
                body=json.dumps({}),
                properties=pika.BasicProperties(
                    delivery_mode=2,  # make message persistent
                ),
            )
            print(
                " [x] Sent message to {}".format(routing_key_done),
                flush=True,
            )
        print(" [x] Done", flush=True)
        ch.basic_ack(delivery_tag=method.delivery_tag)

    channel.basic_qos(prefetch_count=1)
    channel.basic_consume(
        queue="preprocessing_python_task_queue", on_message_callback=callback
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
