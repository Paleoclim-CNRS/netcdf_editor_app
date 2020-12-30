from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, current_app, session
)
from werkzeug.exceptions import abort
from werkzeug.utils import secure_filename

import os
import tempfile
import functools

from netcdf_editor_app.auth import login_required
from netcdf_editor_app.db import get_db

import xarray as xr
import numpy as np

bp = Blueprint('app', __name__)


@bp.route('/')
@login_required
def index():
    # Remove the datafile if we go back to datafile selection screen
    session.pop('data_file_id', None)
    db = get_db()
    data_files = db.execute(
        'SELECT created, filename, username, owner_id, df.id'
        ' FROM data_files df JOIN user u ON df.owner_id = u.id'
        ' ORDER BY created DESC'
    ).fetchall()

    return render_template('app/index.html', data_files=data_files)


def allowed_file(filename):
    ALLOWED_EXTENSIONS = {'nc'}
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@bp.route('/upload', methods=('GET', 'POST'))
@login_required
def upload():
    if request.method == 'POST':
        # check if the post request has the file part
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        # if user does not select file, browser also
        # submit an empty part without filename
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            temp_name = next(tempfile._get_candidate_names()) + ".nc"
            file.save(os.path.join(
                current_app.config['UPLOAD_FOLDER'], temp_name))
            db = get_db()
            data_file = db.execute(
                'INSERT INTO data_files (owner_id, filename, filepath)'
                ' VALUES (?, ?, ?)',
                (g.user['id'], filename, temp_name)
            )
            db.commit()
            data_file_id = data_file.lastrowid
            return redirect(url_for('set_coords', data_file_id=data_file_id))

    return render_template('app/upload.html')

@bp.route('/set_coords', methods=('GET', 'POST'))
@login_required
def set_coords():
    _id = request.args['data_file_id']
    db = get_db()
    if request.method == 'POST':
        lat = request.form['Latitude']
        lon = request.form['Longitude']
        db.execute(
            'UPDATE data_files SET longitude = ?, latitude = ? WHERE id = ?', (lon, lat, _id)
        )
        db.commit()
        return redirect(url_for('index'))
          

    filepath = db.execute(
        'SELECT filepath FROM data_files WHERE id = ?', (request.args['data_file_id'])
    ).fetchone()['filepath']
    filepath = os.path.join(current_app.instance_path, filepath)
    coordinate_names = [name for name in xr.open_dataset(filepath).coords]
    return render_template('app/set_coords.html', coordinate_names = coordinate_names)


@bp.route('/steps')
@login_required
def steps():
    data_file_id_req = request.args.get('id')
    if data_file_id_req is not None:
        session['data_file_id'] = data_file_id_req
    
    if 'data_file_id' not in session.keys() or session['data_file_id'] is None:
        return redirect(url_for('index'))
    data_file_id = session['data_file_id']

    db = get_db()
    data_file_name = db.execute(
        'SELECT filename FROM data_files WHERE id = ?', (data_file_id, )
    ).fetchone()['filename']
    return render_template('app/steps.html', data_file_name = data_file_name )


def data_file_required(view):
    @functools.wraps(view)
    def wrapped_view(**kwargs):
        if 'data_file_id' not in session.keys() or session['data_file_id'] is None:
            return redirect(url_for('index'))

        return view(**kwargs)

    return wrapped_view


@bp.route('/regrid', methods=('GET', 'POST'))
@login_required
@data_file_required
def regrid():
    if request.method == 'POST':
        lon_step = float(request.form['Longitude Step'])
        lat_step = float(request.form['Latitude Step'])
        interpolator = request.form['interpolator']
        error = ''

        if not lon_step or lon_step < 0:
            error += 'Incorrect Longitude step; '
        elif not lat_step or lat_step < 0:
            error += 'Incorrect Latitude step; '
        elif interpolator not in ['linear', 'neareast']:
            error += "Unknown interpolator"

        if not len(error):
            # Get filename
            db = get_db()
            filepath = db.execute(
                'SELECT filepath FROM data_files WHERE id = ?', (session['data_file_id'], )
            ).fetchone()['filepath']
            # Load file
            full_filepath = os.path.join(current_app.instance_path, filepath)
            ds = xr.open_dataset(full_filepath)
            # Interpolate data file
            new_lon = np.arange(ds.lon[0], ds.lon[-1], lon_step)
            new_lat = np.arange(ds.lat[0], ds.lat[-1], lat_step)
            ds = ds.interp(lat=lat_step, lon=lon_step, method=interpolator)
            # Save file
            # return redirect(url_for('app.steps'))
        
        flash(error)

    return render_template('app/regrid.html')
