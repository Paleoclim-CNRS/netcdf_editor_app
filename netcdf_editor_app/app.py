from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, current_app, session
)
from werkzeug.exceptions import abort
from werkzeug.utils import secure_filename

import os
import tempfile
import functools

from netcdf_editor_app.auth import login_required
from netcdf_editor_app.db import load_file, get_file_path, get_lon_lat_names, get_db

import xarray as xr
import numpy as np
import hvplot.xarray

import holoviews as hv
from bokeh.resources import CDN
from bokeh.embed import file_html, server_session, components
from bokeh.embed import file_html

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
            return redirect(url_for('app.set_coords', _id=data_file_id))

    return render_template('app/upload.html')


@bp.route('/<int:_id>/set_coords', methods=('GET', 'POST'))
@login_required
def set_coords(_id):
    db = get_db()
    if request.method == 'POST':
        lat = request.form['Latitude']
        lon = request.form['Longitude']
        db.execute(
            'UPDATE data_files SET longitude = ?, latitude = ? WHERE id = ?', (lon, lat, str(
                _id))
        )
        db.commit()
        return redirect(request.form['next'])

    filepath = db.execute(
        'SELECT filepath FROM data_files WHERE id = ?', (str(_id))
    ).fetchone()['filepath']
    filepath = os.path.join(current_app.instance_path, filepath)
    coordinate_names = [name for name in xr.open_dataset(filepath).coords]
    return render_template('app/set_coords.html', coordinate_names=coordinate_names)


@bp.route('/<int:_id>/steps')
@login_required
def steps(_id):
    db = get_db()
    data_file_name = db.execute(
        'SELECT filename FROM data_files WHERE id = ?', (str(_id), )
    ).fetchone()['filename']
    return render_template('app/steps.html', data_file_name=data_file_name, _id=_id)


@bp.route('/<int:_id>/map')
@login_required
def map(_id):
    ds = load_file(_id)
    lon, lat = get_lon_lat_names(_id)
    plot = ds.hvplot(x=lon, y=lat).opts(responsive=True)
    plot = hv.render(plot, backend='bokeh')
    plot.sizing_mode = 'scale_width'
    script, div = components(plot)
    return render_template('app/map.html', script=script, div=div, data_file_id=_id)


@bp.route('/<int:_id>/regrid', methods=('GET', 'POST'))
@login_required
def regrid(_id):
    if request.method == 'POST':
        lon_step = float(request.form['Longitude Step'])
        lat_step = float(request.form['Latitude Step'])
        interpolator = request.form['interpolator']
        error = ''

        if not lon_step or lon_step < 0:
            error += 'Incorrect Longitude step; '
        elif not lat_step or lat_step < 0:
            error += 'Incorrect Latitude step; '
        elif interpolator not in ['linear', 'nearest']:
            error += "Unknown interpolator"

        if not len(error):
            # Load file
            ds = load_file(_id)
            lon, lat = get_lon_lat_names(_id)
            # Interpolate data file
            new_lon = np.arange(ds[lon][0], ds[lon][-1], lon_step)
            new_lat = np.arange(ds[lat][0], ds[lat][-1], lat_step)
            interp_options = {
                lon: new_lon,
                lat: new_lat,
            }
            ds = ds.interp(interp_options, method=interpolator,)
            # Save file
            ds.to_netcdf(get_file_path(_id))
            flash("File regrided using {} interpolation with Longitude steps {} and Latitude steps {}".format(
                interpolator, lon_step, lat_step))
            return redirect(url_for('app.steps', _id=_id))

        flash(error)

    return render_template('app/regrid.html')


@bp.route('/<int:_id>/internal_oceans')
@login_required
def internal_oceans(_id):
    script = server_document(url='http://localhost:5006/value_changer',
                             arguments={'id': _id})
    # Arguments are reached through Bokeh curdoc.session_context.request.arguments
    # And hence through panel.state.curdoc.session_context.request.arguments
    return render_template("app/internal_oceans.html", script=script)

