from flask import (
    Blueprint,
    flash,
    redirect,
    render_template,
    request,
    url_for,
    session,
    current_app,
    send_from_directory,
)
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
from netcdf_editor_app.utils.pft import generate_pft_netcdf
from netcdf_editor_app.utils.heatflow import create_heatflow
from netcdf_editor_app.utils.ahmcoef import create_ahmcoef

import numpy
import xarray as xr
import pandas as pd
import json
import hvplot.xarray

import holoviews as hv
from bokeh.embed import components, server_document, file_html
from bokeh.resources import CDN
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, CustomJS, Slider, Select
from bokeh.plotting import Figure, output_file, show

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


@bp.route("/<int:_id>/<string:file_type>/download", methods=["GET"])
@login_required
def download(_id, file_type):
    data_file_name = get_filename(_id)
    name, extension = data_file_name.split(".")
    name += "_" + file_type + "_netcdf_flask_app"
    data_file_name = name + "." + extension

    filename = get_file_path(_id, file_type, full=False)
    uploads = os.path.join(current_app.root_path, current_app.config["UPLOAD_FOLDER"])
    return send_from_directory(
        directory=uploads,
        filename=filename,
        as_attachment=True,
        attachment_filename=data_file_name,
    )


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


@bp.route("/<int:_id>/<string:file_type>/view", methods=["GET"])
@login_required
def view_database_file(_id, file_type):

    ds = load_file(_id, file_type)
    d = {}
    for data_var in ds.data_vars:
        d[data_var] = [ds[data_var].values]
    d["to_plot"] = [ds[list(ds.data_vars)[0]].values]
    source = ColumnDataSource(d)

    callback = CustomJS(
        args=dict(source=source),
        code="""
        var data = source.data;
        data['to_plot'] = data[cb_obj.value];
        source.change.emit();
    """,
    )
    select = Select(title="Variable:", options=list(ds.data_vars))
    select.js_on_change("value", callback)

    p = Figure(x_range=(-180, 180), y_range=(-90, 90), aspect_ratio=2.5, tools='pan,wheel_zoom,box_zoom,reset, hover')
    p.sizing_mode = "scale_width"
    p.image(
        image="to_plot",
        x=-180,
        y=-90,
        dw=360,
        dh=180,
        source=source,
        palette="Viridis11",
    )
    script, div = components(column(column(select), p, sizing_mode="stretch_width"))

    return render_template(
        "app/bokeh_plot.html",
        script=script,
        div=div,
        data_file_id=_id,
        title=file_type.capitalize(),
    )


@bp.route("/<int:_id>/steps")
@login_required
def steps(_id):
    data_file_name = get_filename(_id)

    seen_file_types = get_file_types(_id)
    data = []
    for name in seen_file_types:
        data.append(
            [
                name.capitalize(),
                f"<form action=\"{ url_for('app.view_database_file', _id=_id, file_type=name.lower()) }\" method=\"GET\"> \
                    <button type=\"submit\" class=\"btn btn-info\"><i class=\"fas fa-map\"></i> View</button> \
                </form>",
                f"<form action=\"{ url_for('app.variable_explorer', _id=_id, file_type=name.lower()) }\" method=\"GET\"> \
                    <button type=\"submit\" class=\"btn btn-info\"><i class=\"fas fa-columns\"></i> View</button> \
                </form>",
                f"<form action=\"{ url_for('app.revision_comparison', _id=_id, file_type=name.lower()) }\" method=\"GET\"> \
                    <button type=\"submit\" class=\"btn btn-info\"><i class=\"fas fa-arrows-alt-h\"></i> View</button> \
                </form>",
                f"<form action=\"{ url_for('app.download', _id=_id, file_type=name.lower()) }\" method=\"GET\"> \
                    <button type=\"submit\" class=\"btn btn-primary\"><i class=\"fas fa-download\"></i> Download</button> \
                </form>",
            ]
        )

    df = pd.DataFrame(
        data,
        columns=[
            "File Type",
            "View",
            "Complex Viewer",
            "Revision Comparison",
            "Download Link",
        ],
    )
    return render_template(
        "app/steps.html",
        data_file_name=data_file_name,
        _id=_id,
        assets_table=df.to_html(
            index=False, justify="center", border=0, classes="table", escape=False
        ),
    )


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


@bp.route("/<int:_id>/<string:file_type>/revision_comparison")
@login_required
def revision_comparison(_id, file_type):
    try:
        ds_latest = load_file(_id, file_type, -1)
        ds_previous = load_file(_id, file_type, -2)
    except IndexError:
        flash("Not enough revisions")
        return redirect(url_for("app.steps", _id=_id))
    ds = ds_latest - ds_previous
    lon, lat = get_lon_lat_names(_id)
    plot = ds.hvplot(x=lon, y=lat).opts(responsive=True, cmap="terrain")
    plot = hv.render(plot, backend="bokeh")
    plot.sizing_mode = "scale_width"
    script, div = components(plot)
    return render_template("app/map.html", script=script, div=div, data_file_id=_id)


@bp.route("/<int:_id>/<string:file_type>/variable_explorer")
@login_required
def variable_explorer(_id, file_type):
    script = server_document(
        url=f"{url_for('index')}panel/value_changer",
        arguments={
            "id": _id,
            "redirect": url_for("app.steps", _id=_id),
            "file_type": file_type,
        },
        # resources=CDN
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
        limits = request.form['limits']
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
        elif limits not in ["default", "data"]:
            error += "Unknown limit type"

        if not len(error):
            # Load file
            ds = load_file(_id, "raw")
            lon, lat = get_lon_lat_names(_id)
            # Extremities
            new_values = []
            # Limits
            default_limits = [180, 90]
            for coord, step, default_limit  in zip([lon, lat], [lon_step, lat_step], default_limits):
                if limits == 'default':
                    lower = -default_limit
                    upper = default_limit
                elif limits == 'data':
                    # Get sorted values
                    sorted_vals = numpy.sort(numpy.unique(ds[coord]))
                    lower = ds[coord].min() - (sorted_vals[1] - sorted_vals[0]) / 2.0
                    upper = ds[coord].max() + (sorted_vals[-1] - sorted_vals[-2]) / 2.0
                else:
                    raise AttributeError("Unknown data type passed from UI")

                min_val = lower + step / 2.0
                max_val = upper + step / 2.0

                # TODO maybe we should use numpy.linspace here?
                new_values.append(numpy.arange(min_val, max_val, step))
            # Interpolate data file
            interp_options = {
                lon: new_values[0],
                lat: new_values[1],
            }
            ds = ds.interp(
                interp_options,
                method=interpolator,
                kwargs=dict(fill_value=None)
            )
            # Save file
            save_revision(_id, ds, "raw")
            flash(
                "File regrided using {} interpolation with Longitude steps {} and Latitude steps {} and {} limits".format(
                    interpolator, lon_step, lat_step, limits
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
        url=f"{url_for('index')}panel/internal_oceans",
        arguments={"id": _id, "redirect": url_for("app.steps", _id=_id)},
    )
    # Arguments are reached through Bokeh curdoc.session_context.request.arguments
    # And hence through panel.state.curdoc.session_context.request.arguments
    return render_template("app/panel_app.html", script=script, title="Internal Oceans")


@bp.route("/<int:_id>/routing", methods=("GET", "POST"))
@login_required
def routing(_id):
    ds = load_file(_id, "raw")
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
            ds_routing, ds_bathy, ds_soils, ds_topo_high_res = run_routines(
                topography, latitudes
            )
            save_revision(_id, ds_routing, "routing")
            save_revision(_id, ds_bathy, "bathy")
            save_revision(_id, ds_soils, "soils")
            save_revision(_id, ds_topo_high_res, "topo_high_res")

            flash("Routing run succesfully")

            return redirect(url_for("app.steps", _id=_id))

        flash(error)

    data_shape = tuple(ds.dims.values())
    show_regrid = False
    if data_shape != (180, 360):
        show_regrid = True

    return render_template(
        "app/routing.html",
        _id=_id,
        variable_names=variable_names,
        show_regrid=show_regrid,
        data_shape=data_shape,
    )


@bp.route("/<int:_id>/passage_problems")
@login_required
def passage_problems(_id):
    seen_file_types = get_file_types(_id)
    if "routing" not in seen_file_types:
        script = f"<div class=\"alert alert-danger\" role=\"alert\" \
            style=\"display: flex; align-items: center; justify-content: center; flex-direction: column;\">\
            Routing file not found in Database please perform Routing step first\
            <br></br> \
            <button onclick=\"window.location.href = \
                '{url_for('app.routing', _id=_id)}';\" id=\"myButton\" class=\"btn btn-primary\" >Routing</button> \
            </div>"
    else:
        script = server_document(
            url=f"{url_for('index')}panel/passage_problems",
            arguments={"id": _id, "redirect": url_for("app.steps", _id=_id)},
        )
    # Arguments are reached through Bokeh curdoc.session_context.request.arguments
    # And hence through panel.state.curdoc.session_context.request.arguments
    return render_template(
        "app/panel_app.html", script=script, title="Passage Problems"
    )


@bp.route("/<int:_id>/pft", methods=("GET", "POST"))
@login_required
def pft(_id):
    if request.method == "POST":
        data = json.loads(request.form["data"])
        resp_array = numpy.array(data["dataArray"])
        # Make sure the data array is in the expected format
        # First col = cutoff latitudes
        # Next 13 cols are pft types
        assert len(resp_array[0] == 14)
        pft_values = resp_array[:, 1:]
        latitudes = resp_array[:, 0]
        # Make sure 90 is the last value
        assert latitudes[-1] == 90

        # Load routing file with final topography
        ds = load_file(_id, "routing")
        assert set(ds.dims) == set(("x", "y"))
        assert len(ds.coords) == 2
        # The PFT values are on a 360 x 720 grid
        # So we need to interpolate the values onto this grid
        lat_vals = numpy.arange(0, 180, 0.5)
        lon_vals = numpy.arange(0, 360, 0.5)
        ds = ds.interp({"y": lat_vals, "x": lon_vals})
        topo = ds.topo.values

        ds = generate_pft_netcdf(topo, latitudes, pft_values)
        save_revision(_id, ds, "pft")
        return redirect(url_for("app.steps", _id=_id))

    seen_file_types = get_file_types(_id)
    not_seen = False
    if "routing" not in seen_file_types:
        not_seen=True

    return render_template("app/pft.html", _id=_id, not_seen=not_seen)


@bp.route("/<int:_id>/sub_basins")
@login_required
def subbasins(_id):
    script = server_document(
        url=f"{url_for('index')}panel/sub_basins",
        arguments={"id": _id, "redirect": url_for("app.steps", _id=_id)},

    )
    # Arguments are reached through Bokeh curdoc.session_context.request.arguments
    # And hence through panel.state.curdoc.session_context.request.arguments
    return render_template("app/panel_app.html", script=script, title="Sub Basins")

@bp.route("/<int:_id>/heatflow", methods=("GET", "POST"))
@login_required
def heatflow(_id):
    if request.method == 'POST':
        ds = load_file(_id, "bathy")
        
        ds_out = create_heatflow(ds)
        save_revision(_id, ds_out, 'heatflow')

        ds_out = create_ahmcoef(ds)
        save_revision(_id, ds_out, 'ahmcoef')

        return redirect(url_for("app.steps", _id=_id))
        
        flash(error)
    show_routing = "bathy" not in get_file_types(_id)
    return render_template("app/heatflow.html", _id=_id, show_routing=show_routing)