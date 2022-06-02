from datetime import datetime
import numpy
import json

from climate_simulation_platform import create_app
from climate_simulation_platform.db import (
    get_lon_lat_names,
    load_file,
    save_revision,
    step_seen,
    invalidate_step,
)
from climpy.bc.ipsl.routing import run_routines
from climpy.bc.ipsl.heatflow import create_heatflow
from climpy.bc.ipsl.ahmcoef import create_ahmcoef
from climpy.bc.ipsl.pft import generate_pft_netcdf


def regrid(body):
    app = create_app()

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
    print(f" [x] {datetime.now()} interpolation finished", flush=True)
    # Save file
    with app.app_context():
        save_revision(_id, ds, "raw")


def routing(body):
    app = create_app()

    topo_variable = body["topo_var"]
    _id = body["id"]

    save_routing = True
    save_bathy = True
    save_soils = True
    save_topo_high_res = True

    if body.get("next_step_only", False):
        save_bathy = False

    # Load file
    with app.app_context():
        ds = load_file(_id, "raw")
        ds_orca = load_file(_id, "coords")
        lon, lat = get_lon_lat_names(_id)

    latitudes = ds[lat].values
    topography = ds[topo_variable].values
    ds_routing, ds_bathy, ds_soils, ds_topo_high_res = run_routines(
        topography, latitudes, ds_orca
    )
    with app.app_context():
        if save_routing:
            save_revision(_id, ds_routing, "routing")
        if save_bathy:
            save_revision(_id, ds_bathy, "bathy")
        if save_soils:
            save_revision(_id, ds_soils, "soils")
        if save_topo_high_res:
            save_revision(_id, ds_topo_high_res, "topo_high_res")


def pft(body):
    app = create_app()
    _id = body["id"]

    data = json.loads(body["data"])
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
    with app.app_context():
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
    with app.app_context():
        save_revision(_id, ds, "pft")


def heatflow(body):
    app = create_app()
    _id = body["id"]

    with app.app_context():
        ds = load_file(_id, "bathy")

    ds_out = create_heatflow(ds)

    with app.app_context():
        save_revision(_id, ds_out, "heatflow")


def ahmcoef(body):
    app = create_app()
    _id = body["id"]

    with app.app_context():
        ds = load_file(_id, "bathy")

    ds_out = create_ahmcoef(ds)

    with app.app_context():
        save_revision(_id, ds_out, "ahmcoef")


def invalidate(body):
    app = create_app()
    _id = body["id"]

    with app.app_context():
        for task in body["tasks"]:
            if step_seen(_id, task):
                invalidate_step(_id, task)
