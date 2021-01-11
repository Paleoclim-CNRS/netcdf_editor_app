from flask import (
    Blueprint, flash, redirect, render_template, request, url_for, session
)


from netcdf_editor_app.auth import login_required
from netcdf_editor_app.db import get_coord_names, get_db, get_latest_file_versions, get_lon_lat_names, load_file, remove_data_file, save_revision, set_data_file_coords, upload_file

import numpy as np
import hvplot.xarray

import holoviews as hv
from bokeh.embed import components, server_document

bp = Blueprint('app', __name__)


@bp.route('/')
@login_required
def index():
    # Remove the datafile if we go back to datafile selection screen
    session.pop('data_file_id', None)
    data_files = get_latest_file_versions()

    return render_template('app/datafiles.html', data_files=data_files)


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
            data_file_id = upload_file(file)
            return redirect(url_for('app.set_coords', _id=data_file_id))

    return render_template('app/upload.html')


@bp.route('/<int:_id>/delete', methods=('GET', 'POST'))
@login_required
def delete(_id):
    if request.method == 'POST':
        remove_data_file(_id)
        flash("File deleted with id: {}".format(_id))
    return redirect(url_for('index'))


@bp.route('/<int:_id>/set_coords', methods=('GET', 'POST'))
@login_required
def set_coords(_id):
    if request.method == 'POST':
        lat = request.form['Latitude']
        lon = request.form['Longitude']
        set_data_file_coords(_id, longitude=lon, latitude=lat)

        next_page = request.form['next']
        if '/upload' in next_page:
            next_page = url_for('index')
        return redirect(next_page)

    coordinate_names = get_coord_names(_id)
    return render_template('app/set_coords.html', coordinate_names=coordinate_names)


@bp.route('/<int:_id>/steps')
@login_required
def steps(_id):
    db = get_db()
    data_file_name = db.execute(
        'SELECT filename FROM data_files WHERE id = ?', (str(_id), )
    ).fetchone()['filename']
    return render_template('app/steps.html', data_file_name=data_file_name, _id=_id)


@bp.route('/<int:_id>')
@login_required
def redirect_steps(_id):
    return redirect(url_for('app.steps', _id=_id))


@bp.route('/<int:_id>/map')
@login_required
def map(_id):
    ds = load_file(_id)
    lon, lat = get_lon_lat_names(_id)
    plot = ds.hvplot(x=lon, y=lat).opts(responsive=True, cmap='terrain')
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
            # Extremities
            new_values = []
            for coord, step in zip([lon, lat], [lon_step, lat_step]):
                # Get sorted values
                sorted_vals = np.sort(np.unique(ds[coord]))
                min_val = ds[coord].min() - (sorted_vals[1] -
                                             sorted_vals[0]) / 2. + step / 2.
                max_val = ds[coord].max() + (sorted_vals[-1] -
                                             sorted_vals[-2]) / 2. + step / 2.
                new_values.append(np.arange(min_val, max_val, step))
            # Interpolate data file
            interp_options = {
                lon: new_values[0],
                lat: new_values[1],
            }
            ds = ds.interp(interp_options, method=interpolator,)
            # Save file
            save_revision(_id, ds)
            flash("File regrided using {} interpolation with Longitude steps {} and Latitude steps {}".format(
                interpolator, lon_step, lat_step))
            return redirect(url_for('app.steps', _id=_id))

        flash(error)

    return render_template('app/regrid.html')


@bp.route('/<int:_id>/internal_oceans')
@login_required
def internal_oceans(_id):
    script = server_document(url='http://localhost:5006/internal_oceans',
                             arguments={'id': _id})
    # Arguments are reached through Bokeh curdoc.session_context.request.arguments
    # And hence through panel.state.curdoc.session_context.request.arguments
    return render_template("app/panel_app.html", script=script, title="Internal Oceans")


@bp.route('/<int:_id>/passage_problems')
@login_required
def passage_problems(_id):
    script = server_document(url='http://localhost:5006/passage_problems',
                             arguments={'id': _id})
    # Arguments are reached through Bokeh curdoc.session_context.request.arguments
    # And hence through panel.state.curdoc.session_context.request.arguments
    return render_template("app/panel_app.html", script=script, title="Passage Problems")
