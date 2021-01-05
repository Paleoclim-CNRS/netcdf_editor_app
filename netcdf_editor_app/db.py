import sqlite3

import click
from flask import current_app, g
from flask.cli import with_appcontext

import xarray as xr
import os


def get_db():
    if 'db' not in g:
        g.db = sqlite3.connect(
            current_app.config['DATABASE'],
            detect_types=sqlite3.PARSE_DECLTYPES
        )
        g.db.row_factory = sqlite3.Row

    return g.db


def close_db(e=None):
    db = g.pop('db', None)

    if db is not None:
        db.close()


def init_db():
    db = get_db()

    with current_app.open_resource('db_schema.sql') as f:
        db.executescript(f.read().decode('utf8'))


def load_file(_id):
    # Get filename
    filepath = get_file_path(_id)
    # Load file
    return xr.open_dataset(filepath)


def get_file_path(_id, full=True):
    db = get_db()
    filepath = db.execute(
        'SELECT filepath FROM data_files WHERE id = ?', (str(_id), )
    ).fetchone()['filepath']
    if not full:
        return filepath
    return os.path.join(current_app.instance_path, filepath)


def get_lon_lat_names(_id):
    db = get_db()
    return db.execute(
        'SELECT longitude, latitude FROM data_files WHERE id = ?', (str(_id), )
    ).fetchone()


@click.command('init-db')
@with_appcontext
def init_db_command():
    """Clear the existing data and create new tables."""
    init_db()
    click.echo('Initialized the database.')


def init_app(app):
    app.teardown_appcontext(close_db)
    app.cli.add_command(init_db_command)
