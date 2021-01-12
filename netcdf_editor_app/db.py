import sqlite3
import os
import tempfile

import click
from flask import current_app, g
from flask.cli import with_appcontext

from werkzeug.utils import secure_filename

import xarray as xr


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


def load_file(_id, revision=-1):
    # Get filename
    filepath = get_file_path(_id, revision=revision)
    # Load file
    ds = xr.open_dataset(filepath)
    return ds


def get_coord_names(_id):
    ds = load_file(_id)
    return [name for name in ds.coords]


def set_data_file_coords(_id, longitude, latitude):
    db = get_db()
    db.execute(
        'UPDATE data_files SET longitude = ?, latitude = ? WHERE id = ?', (longitude, latitude, str(
            _id))
    )
    db.commit()


def upload_file(file):
    filename = secure_filename(file.filename)
    temp_name = next(tempfile._get_candidate_names()) + ".nc"
    # Save the file to the file system
    file.save(os.path.join(
        current_app.config['UPLOAD_FOLDER'], temp_name))
    # Add file to known data_files
    db = get_db()
    data_file = db.execute(
        'INSERT INTO data_files (owner_id, filename)'
        ' VALUES (?, ?)',
        (g.user['id'], filename)
    )
    db.commit()
    data_file_id = data_file.lastrowid
    # ADD the file to the revisions table
    db.execute(
        'INSERT INTO revisions (data_file_id, filepath, revision)'
        ' VALUES (?, ?, ?)',
        (data_file_id, temp_name, 0)
    )
    db.commit()

    return data_file_id


def save_revision(_id, ds):
    temp_name = next(tempfile._get_candidate_names()) + ".nc"
    # Save the file to the file system
    ds.to_netcdf(os.path.join(
        current_app.config['UPLOAD_FOLDER'], temp_name))
    # ADD the file to the revisions table
    db = get_db()
    # Get latest revision
    query = 'SELECT revision FROM revisions WHERE data_file_id = ? ORDER BY revision DESC LIMIT 0, 1'
    revision_nb = db.execute(query, (str(_id),)).fetchone()['revision'] + 1
    db.execute(
        'INSERT INTO revisions (data_file_id, filepath, revision)'
        ' VALUES (?, ?, ?)',
        (str(_id), temp_name, revision_nb)
    )
    db.commit()


def get_latest_file_versions():
    query = 'SELECT created, filename, df.id FROM' +\
        ' (SELECT MAX(revision) as revision, created, data_file_id FROM revisions GROUP BY data_file_id) as r' +\
        ' JOIN data_files df ON r.data_file_id = df.id' +\
        ' JOIN user u ON df.owner_id = u.id' +\
        ' ORDER BY created DESC'
    db = get_db()
    data_files = db.execute(query).fetchall()
    return data_files


def get_file_path(_id, full=True, revision=-1):
    db = get_db()
    revisions = db.execute(
        'SELECT revision FROM revisions WHERE data_file_id = ? ORDER BY revision ASC', (str(
            _id), )
    ).fetchall()
    revisions = [rev['revision'] for rev in revisions]
    revision_nb = revisions[revision]
    filepath = db.execute(
        'SELECT filepath FROM revisions WHERE data_file_id = ? AND revision = ?', (str(
            _id), str(revision_nb), )
    ).fetchone()['filepath']
    if not full:
        return filepath
    return os.path.join(current_app.instance_path, filepath)


def get_all_file_paths(_id, full=True):
    db = get_db()
    filepath = db.execute(
        'SELECT filepath FROM revisions WHERE data_file_id = ?', (str(_id), )
    ).fetchall()
    filepath = [fp['filepath'] for fp in filepath]
    if not full:
        return filepath
    return [os.path.join(current_app.instance_path, fp) for fp in filepath]


def remove_data_file(_id):
    # Get filepath
    file_paths = get_all_file_paths(_id)
    for filepath in file_paths:
        os.remove(filepath)
    # Remove entry from DB
    db = get_db()
    db.execute(
        'DELETE FROM data_files WHERE id = ?', (str(_id), )
    )
    db.execute(
        'DELETE FROM revisions WHERE data_file_id = ?', (str(_id), )
    )
    db.commit()


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
