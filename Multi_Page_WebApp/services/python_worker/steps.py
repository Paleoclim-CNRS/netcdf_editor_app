import numpy

from netcdf_editor_app import create_app
from netcdf_editor_app.db import get_lon_lat_names, load_file, save_revision


def regrid(body):
    app = create_app()
    print(" [x] Body in function %r" % body, flush=True)

    limits = body["limits"]
    lon_step = float(body["Longitude Step"])
    lat_step = float(body["Latitude Step"])
    interpolator = body["interpolator"]
    _id = body["id"]

    # Load file
    with app.app_context():
        ds = load_file(_id, "raw")
        lon, lat = get_lon_lat_names(_id)
    # Extremities
    new_values = []
    # Limits
    default_limits = [180, 90]
    for coord, step, default_limit in zip(
        [lon, lat], [lon_step, lat_step], default_limits
    ):
        if limits == "default":
            lower = -default_limit
            upper = default_limit
        elif limits == "data":
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
    ds = ds.interp(interp_options, method=interpolator, kwargs=dict(fill_value=None))
    print(" [x] interpolation finished", flush=True)
    # Save file
    with app.app_context():
        save_revision(_id, ds, "raw")


def routing(body):
    app = create_app()

    topo_variable = body["topo_var"]
    _id = body["id"]

    # Load file
    with app.app_context():
        ds = load_file(_id, "raw")
        lon, lat = get_lon_lat_names(_id)

    latitudes = ds[lat].values
    topography = ds[topo_variable].values
    ds_routing, ds_bathy, ds_soils, ds_topo_high_res = run_routines(
        topography, latitudes
    )
    with app.app_context():
        save_revision(_id, ds_routing, "routing")
        save_revision(_id, ds_bathy, "bathy")
        save_revision(_id, ds_soils, "soils")
        save_revision(_id, ds_topo_high_res, "topo_high_res")
