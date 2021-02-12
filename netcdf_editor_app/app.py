from flask import Blueprint, flash, redirect, render_template, request, url_for, session, current_app, send_from_directory
import werkzeug

import os

from netcdf_editor_app.auth import login_required
from netcdf_editor_app.db import (
    get_coord_names,
    get_latest_file_versions,
    get_lon_lat_names,
    get_file_path,
    get_file_types,
    load_file,
    remove_data_file,
    save_revision,
    set_data_file_coords,
    upload_file,
    get_filename,
)

from netcdf_editor_app.utils.routing import run_routines

import numpy
import pandas as pd
import hvplot.xarray

import holoviews as hv
from bokeh.embed import components, server_document

bp = Blueprint("app", __name__)


@bp.route("/")
@login_required
def index():
    # Remove the datafile if we go back to datafile selection screen
    session.pop("data_file_id", None)
    data_files = get_latest_file_versions()

    return render_template("app/datafiles.html", data_files=data_files)


def allowed_file(filename):
    ALLOWED_EXTENSIONS = {"nc"}
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS


@bp.route("/upload", methods=("GET", "POST"))
@login_required
def upload():
    if request.method == "POST":
        # check if the post request has the file part
        if "file" not in request.files:
            flash("No file part")
            return redirect(request.url)
        file = request.files["file"]
        # if user does not select file, browser also
        # submit an empty part without filename
        if file.filename == "":
            flash("No selected file")
            return redirect(request.url)
        if file and allowed_file(file.filename):
            data_file_id = upload_file(file)
            return redirect(url_for("app.set_coords", _id=data_file_id))

    return render_template("app/upload.html")


@bp.route('/<int:_id>/<string:file_type>/download', methods=['GET'])
@login_required
def download(_id, file_type):
    data_file_name = get_filename(_id)
    name, extension = data_file_name.split(".")
    name += "_" + file_type + "_netcdf_flask_app"
    data_file_name = name + "." + extension

    filename = get_file_path(_id, file_type, full=False)
    uploads = os.path.join(current_app.root_path, current_app.config['UPLOAD_FOLDER'])
    return send_from_directory(directory=uploads, filename=filename, as_attachment=True, attachment_filename=data_file_name)


@bp.route("/<int:_id>/delete", methods=("GET", "POST"))
@login_required
def delete(_id):
    if request.method == "POST":
        remove_data_file(_id)
        flash("File deleted with id: {}".format(_id))
    return redirect(url_for("index"))


@bp.route("/<int:_id>/set_coords", methods=("GET", "POST"))
@login_required
def set_coords(_id):
    if request.method == "POST":
        lat = request.form["Latitude"]
        lon = request.form["Longitude"]
        set_data_file_coords(_id, longitude=lon, latitude=lat)

        next_page = request.form["next"]
        if "/upload" in next_page:
            next_page = url_for("index")
        return redirect(next_page)

    coordinate_names = get_coord_names(_id)
    return render_template("app/set_coords.html", coordinate_names=coordinate_names)


@bp.route("/<int:_id>/steps")
@login_required
def steps(_id):
    data_file_name = get_filename(_id)

    seen_file_types = get_file_types(_id)
    data = []
    for name in seen_file_types:
        data.append(
            [name.capitalize(),
            f"<form action=\"{ url_for('app.download', _id=_id, file_type=name.lower()) }\" method=\"GET\"> \
                <button type=\"submit\" class=\"btn btn-primary\"><i class=\"fas fa-download\"></i> Download</button> \
            </form>"
            ]
        )

    df = pd.DataFrame(data, columns=["File Type", "Download Link"])
    return render_template("app/steps.html", data_file_name=data_file_name, _id=_id, assets_table=df.to_html(index=False, justify='center', border=0, classes="table", escape=False))


@bp.route("/<int:_id>")
@login_required
def redirect_steps(_id):
    return redirect(url_for("app.steps", _id=_id))


@bp.route("/<int:_id>/map")
@login_required
def map(_id):
    ds = load_file(_id)
    lon, lat = get_lon_lat_names(_id)
    plot = ds.hvplot(x=lon, y=lat).opts(responsive=True, cmap="terrain")
    plot = hv.render(plot, backend="bokeh")
    plot.sizing_mode = "scale_width"
    script, div = components(plot)
    return render_template("app/map.html", script=script, div=div, data_file_id=_id)


@bp.route("/<int:_id>/revision_comparison")
@login_required
def revision_comparison(_id):
    ds_latest = load_file(_id, -1)
    ds_previous = load_file(_id, -2)
    ds = ds_latest - ds_previous
    lon, lat = get_lon_lat_names(_id)
    plot = ds.hvplot(x=lon, y=lat).opts(responsive=True, cmap="terrain")
    plot = hv.render(plot, backend="bokeh")
    plot.sizing_mode = "scale_width"
    script, div = components(plot)
    return render_template("app/map.html", script=script, div=div, data_file_id=_id)


@bp.route("/<int:_id>/variable_explorer")
@login_required
def variable_explorer(_id):
    script = server_document(
        url="http://{os.environ['PANEL_HOST']}:{os.environ['PANEL_SOCKET']}/value_changer", arguments={"id": _id}
    )
    # Arguments are reached through Bokeh curdoc.session_context.request.arguments
    # And hence through panel.state.curdoc.session_context.request.arguments
    return render_template(
        "app/panel_app.html", script=script, title="Variable Explorer"
    )


@bp.route("/<int:_id>/regrid", methods=("GET", "POST"))
@login_required
def regrid(_id):
    if request.method == "POST":
        lon_step = float(request.form["Longitude Step"])
        lat_step = float(request.form["Latitude Step"])
        interpolator = request.form["interpolator"]
        error = ""

        if not lon_step or lon_step < 0:
            error += "Incorrect Longitude step; "
        elif not lat_step or lat_step < 0:
            error += "Incorrect Latitude step; "
        elif interpolator not in ["linear", "nearest"]:
            error += "Unknown interpolator"

        if not len(error):
            # Load file
            ds = load_file(_id, 'raw')
            lon, lat = get_lon_lat_names(_id)
            # Extremities
            new_values = []
            for coord, step in zip([lon, lat], [lon_step, lat_step]):
                # Get sorted values
                sorted_vals = numpy.sort(numpy.unique(ds[coord]))
                min_val = (
                    ds[coord].min()
                    - (sorted_vals[1] - sorted_vals[0]) / 2.0
                    + step / 2.0
                )
                max_val = (
                    ds[coord].max()
                    + (sorted_vals[-1] - sorted_vals[-2]) / 2.0
                    + step / 2.0
                )
                new_values.append(numpy.arange(min_val, max_val, step))
            # Interpolate data file
            interp_options = {
                lon: new_values[0],
                lat: new_values[1],
            }
            ds = ds.interp(
                interp_options,
                method=interpolator,
            )
            # Save file
            save_revision(_id, ds, 'raw')
            flash(
                "File regrided using {} interpolation with Longitude steps {} and Latitude steps {}".format(
                    interpolator, lon_step, lat_step
                )
            )
            try:
                next_page = request.form["next"]
            except werkzeug.exceptions.BadRequestKeyError:
                next_page = url_for("app.steps", _id=_id)
            return redirect(next_page)

        flash(error)

    return render_template("app/regrid.html")


@bp.route("/<int:_id>/internal_oceans")
@login_required
def internal_oceans(_id):
    script = server_document(
        url=f"http://{os.environ['PANEL_HOST']}:{os.environ['PANEL_SOCKET']}/internal_oceans", arguments={"id": _id}
    )
    # Arguments are reached through Bokeh curdoc.session_context.request.arguments
    # And hence through panel.state.curdoc.session_context.request.arguments
    return render_template("app/panel_app.html", script=script, title="Internal Oceans")

@bp.route("/<int:_id>/routing", methods=("GET", "POST"))
@login_required
def routing(_id):
    ds = load_file(_id, 'raw')
    variable_names = list(ds.data_vars)
    lon, lat = get_lon_lat_names(_id)

    if request.method == "POST":
        topo_variable = request.form["topo_var"]
        error = ""

        if not len(topo_variable):
            error += "Topography Variable not understood; "
        elif topo_variable not in variable_names:
            error += "Topography Variable not in data set"

        if not len(error):
            # Load file
            lon, lat = get_lon_lat_names(_id)
            latitudes = ds[lat].values
            topography = ds[topo_variable].values
            ds_routing, ds_bathy, ds_soils, ds_topo_high_res = run_routines(topography, latitudes)
            save_revision(_id, ds_routing, 'routing')
            save_revision(_id, ds_bathy, 'bathy')
            save_revision(_id, ds_soils, 'soils')
            save_revision(_id, ds_topo_high_res, 'topo_high_res')

            flash("Routing run succesfully")

            return redirect(url_for("app.steps", _id=_id))

        flash(error)

    data_shape = tuple(ds.dims.values())
    show_regrid = False
    if data_shape != (180, 360):
        show_regrid = True
    
    return render_template("app/routing.html", _id=_id, variable_names=variable_names, show_regrid=show_regrid, data_shape=data_shape)

@bp.route("/<int:_id>/passage_problems")
@login_required
def passage_problems(_id):
    seen_file_types = get_file_types(_id)
    if "routing" not in seen_file_types:
        script = f"<div class=\"alert alert-danger\" role=\"alert\" \
            style=\"display: flex; align-items: center; justify-content: center; flex-direction: column;\">\
            Routing file not found in Database please perform Routing step first\
            <br></br> \
            <button onclick=\"window.location.href = '{url_for('app.routing', _id=_id)}';\" id=\"myButton\" class=\"btn btn-primary\" >Routing</button> \
            </div>"
    else:
        script = server_document(
            url="http://{os.environ['PANEL_HOST']}:{os.environ['PANEL_SOCKET']}/passage_problems", arguments={"id": _id}
        )
    # Arguments are reached through Bokeh curdoc.session_context.request.arguments
    # And hence through panel.state.curdoc.session_context.request.arguments
    return render_template(
        "app/panel_app.html", script=script, title="Passage Problems"
    )

@bp.route("/<int:_id>/pft",  methods=("GET", "POST"))
@login_required
def pft(_id):
    return render_template("app/pft.html")