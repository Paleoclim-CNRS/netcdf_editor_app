#!/usr/bin/env python
import pika, sys, os
import time
import json
import os

def main():
    connection = pika.BlockingConnection(pika.ConnectionParameters(host=os.environ['BROKER_HOSTNAME']))
    channel = connection.channel()

    channel.queue_declare(queue='task_queue', durable=True)

    def callback(ch, method, properties, body):
        print(" [x] Received %r" % body.decode(), flush=True)
        print(" [x] Received %r" % json.loads(body)['data'], flush=True)
        time.sleep(body.count(b'.'))
        print(" [x] Done", flush=True)
        ch.basic_ack(delivery_tag = method.delivery_tag)

    channel.basic_qos(prefetch_count=1)
    channel.basic_consume(queue='task_queue', on_message_callback=callback)

    print(' [*] Waiting for messages. To exit press CTRL+C', flush=True)
    channel.start_consuming()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('Interrupted')
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
