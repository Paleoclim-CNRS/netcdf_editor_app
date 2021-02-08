from netcdf_editor_app.utils import lmdzor

import numpy
import xarray as xr

import pytest

# This is the input file which holds the base topography
ds_input = xr.open_dataset("./tests/data/input.nc")

# This file has the cmsk to test against and the corrected
# Trip values (rtm)
ds_runoff = xr.open_dataset("./tests/data/runoff.nc")

# Validation for:
# nav_lat
# nav_lon
# orog
# basins
# flength 
ds_pre_stn = xr.open_dataset("./tests/data/pre_stn.nc")

ds_stn = xr.open_dataset("./tests/data/stn.nc")


def nan_equal(a,b):
    try:
        numpy.testing.assert_equal(a,b)
    except AssertionError:
        return False
    return True


def test_trip_values():
    assert numpy.all(lmdzor.trip_values == numpy.array([1, 2, 3, 4, 5, 6, 7, 8]))


def test_next_down_stream_cell():
    assert numpy.all(
        lmdzor.next_downstream_cell
        == numpy.array(
            [
                # 0  1  2  3   4   5   6   7   8
                [0, 1, 1, 0, -1, -1, -1, 0, 1],
                [0, 0, 1, 1, 1, 0, -1, -1, -1],
            ]
        ).T
    )


@pytest.mark.parametrize(
    ("i", "j", "arr", "res"),
    (
        (
            1,
            1,
            numpy.empty((3, 3)),
            numpy.array(
                [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]
            ),
        ),
        (
            0,
            0,
            numpy.empty((3, 3)),
            numpy.array(
                [[0, 2], [0, 0], [0, 1], [0, 2], [0, 0], [0, 1], [1, 2], [1, 0], [1, 1]]
            ),
        ),
        (
            2,
            2,
            numpy.empty((3, 3)),
            numpy.array(
                [[1, 1], [1, 2], [1, 0], [2, 1], [2, 2], [2, 0], [2, 1], [2, 2], [2, 0]]
            ),
        ),
    ),
)
def test_get_adjacent_coords(i, j, arr, res):
    adjacent_coords = lmdzor.get_adjacent_coords(i, j, arr)
    assert numpy.all(adjacent_coords == res)


@pytest.mark.parametrize(
    ("i", "j", "arr", "res"),
    (
        (
            1,
            1,
            numpy.arange(9).reshape(3, 3),
            numpy.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]]),
        ),
        (
            0,
            0,
            numpy.arange(9).reshape(3, 3),
            numpy.array([[2, 0, 1], [2, 0, 1], [5, 3, 4]]),
        ),
        (
            2,
            2,
            numpy.arange(9).reshape(3, 3),
            numpy.array([[4, 5, 3], [7, 8, 6], [7, 8, 6]]),
        ),
    ),
)
def test_get_adjacent_values(i, j, arr, res):
    adjacent_values = lmdzor.get_adjacent_values(i, j, arr)
    assert numpy.all(adjacent_values == res)


@pytest.mark.parametrize(
    ("arr", "ensure_gradient", "res"),
    (
        (
            numpy.arange(4, dtype=float).reshape(2, 2),
            False,
            numpy.array([[1, 0, 1, 0], [1, 0, 1, 0], [3, 2, 3, 2], [3, 2, 3, 2]]),
        ),
        (
            numpy.arange(4, dtype=float).reshape(2, 2),
            True,
            numpy.array([[2, 1, 2, 1], [1, 0, 1, 0], [3, 2, 3, 2], [4, 3, 4, 3]]),
        ),
    ),
)
def test_get_padded_array(arr, ensure_gradient, res):
    padded = lmdzor.get_padded_array(arr, ensure_gradient)
    assert numpy.all(padded == res)


def test_cmsk():
    topo, latitudes = ds_input.topo.values, ds_input.lat.values
    topo = lmdzor.fix_topo(topo, latitudes)
    omsk = lmdzor.calculate_omsk(topo)
    cmsk = lmdzor.calculate_cmsk(omsk)
    assert numpy.all(cmsk == ds_runoff.cmsk.values)

def assert_curvilinear_coordinates(rlon, rlat):
    assert numpy.all(rlon == ds_pre_stn.nav_lon)
    assert numpy.all(rlat == ds_pre_stn.nav_lat)

def assert_basins(basins):
    dicts = []
    for basin_values in [basins, ds_pre_stn.basins.values]:
        vals, counts = numpy.unique(basin_values[numpy.where(~numpy.isnan(basin_values))], return_counts=True)
        d = {}
        for v, c in zip(vals, counts):
            d[v] = c
        dicts.append(d.copy())

    # Make sure there are the same number of basins
    assert numpy.all(dicts[0].keys() == dicts[1].keys())
    # Make sure the size of each basin is the same (the numbers may vary
    # depending on the order they are seen)
    for key in dicts[0].keys():
        assert dicts[0][key] == dicts[1][key]

def assert_river_lengths(river_length):
    assert numpy.nanmax(numpy.abs(river_length - ds_pre_stn.flength.values[::-1])) < 10e-3

def assert_distbox(distbox):
    assert  numpy.nanmax(numpy.abs(distbox - ds_stn.riverl.values[::-1])) < 10e3

def assert_updated_trip(trip):
    assert nan_equal(ds_stn.trip.values[::-1], trip)

def assert_dzz(dzz):
    assert numpy.nanmax(numpy.abs(ds_stn.hdiff.values[::-1] - dzz)) < 100

def assert_topo_index(topo_index):
    assert numpy.nanmax(numpy.abs(ds_stn.topoind.values[::-1] - topo_index)) < 10e4

def test_pre_stn():
    # use trip values provides to test if from trip values we get same values
    trip = ds_runoff.rtm.values
    topo = ds_runoff.topo.values

    omsk = lmdzor.calculate_omsk(topo)

    rlon, rlat = lmdzor.calculate_curvilinear_coordinates()
    assert_curvilinear_coordinates(rlon, rlat)

    area = lmdzor.calculate_area(rlat)
    basins = lmdzor.calculate_basins(topo, trip, area)
    assert_basins(basins)

    ocean_distances = lmdzor.calculate_ocean_distances(trip, rlat)
    river_length = lmdzor.calculate_river_lengths(topo, trip, ocean_distances, omsk)
    assert_river_lengths(river_length)

    distbox = lmdzor.calculate_distbox(river_length, trip)
    assert_distbox(distbox)

    outflow_points = lmdzor.calculate_outflow_points(basins, trip)
    trip = lmdzor.calculate_trip_outflow_values(trip, outflow_points, basins, omsk, rlat)
    assert_updated_trip(trip)

    dzz = lmdzor.calculate_dzz(topo, trip, distbox, omsk)
    assert_dzz(dzz)

    topo_index = lmdzor.calculate_topo_index(distbox, dzz, omsk)
    assert_topo_index(topo_index)
