from netcdf_editor_app.utils import lmdzor

import numpy
import xarray as xr

import pytest

# This is the input file which holds the base topography
ds_input = xr.open_dataset('./tests/data/input.nc')

# This file has the cmsk to test against and the corrected 
# Trip values (rtm)
ds_runoff = xr.open_dataset('./tests/data/runoff.nc')

def test_trip_values():
    assert numpy.all(lmdzor.trip_values == numpy.array([1,2,3,4,5,6,7,8]))
    
def test_next_down_stream_cell():
    assert numpy.all(lmdzor.next_downstream_cell == numpy.array([
    #0  1  2  3   4   5   6   7   8
    [0, 1, 1, 0, -1, -1, -1,  0,  1],
    [0, 0, 1, 1,  1,  0, -1, -1, -1]
]).T)

@pytest.mark.parametrize(('i', 'j', 'arr', 'res'), (
    (1, 1, numpy.empty((3,3)), numpy.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]])),
    (0, 0, numpy.empty((3,3)), numpy.array([[0, 2], [0, 0], [0, 1], [0, 2], [0, 0], [0, 1], [1, 2], [1, 0], [1, 1]])),
    (2, 2, numpy.empty((3,3)), numpy.array([[1, 1], [1, 2], [1, 0], [2, 1], [2, 2], [2, 0], [2, 1], [2, 2], [2, 0]])),
))
def test_get_adjacent_coords(i, j, arr, res):
    adjacent_coords = lmdzor.get_adjacent_coords(i, j, arr)
    assert numpy.all(adjacent_coords == res)

@pytest.mark.parametrize(('i', 'j', 'arr', 'res'), (
    (1, 1, numpy.arange(9).reshape(3,3), numpy.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])), 
    (0, 0, numpy.arange(9).reshape(3,3), numpy.array([[2, 0, 1], [2, 0, 1], [5, 3, 4]])),
    (2, 2, numpy.arange(9).reshape(3,3), numpy.array([[4, 5, 3], [7, 8, 6], [7, 8, 6]]))
))
def test_get_adjacent_values(i, j, arr, res):
    adjacent_values = lmdzor.get_adjacent_values(i, j, arr)
    assert numpy.all(adjacent_values == res)

@pytest.mark.parametrize(('arr', 'ensure_gradient', 'res'), (
    (numpy.arange(4, dtype=float).reshape(2,2), False, numpy.array([[1, 0, 1, 0], [1, 0, 1, 0], [3, 2, 3, 2], [3, 2, 3, 2]])),
    (numpy.arange(4, dtype=float).reshape(2,2), True, numpy.array([[2, 1, 2, 1], [1, 0, 1, 0], [3, 2, 3, 2], [4, 3, 4, 3]])),
))
def test_get_padded_array(arr, ensure_gradient, res):
    padded = lmdzor.get_padded_array(arr, ensure_gradient)
    assert numpy.all(padded == res)

def test_cmsk():
    topo, latitudes = ds_input.topo.values, ds_input.lat.values
    topo = lmdzor.fix_topo(topo, latitudes)
    orog = lmdzor.calculate_orog(topo)
    omsk = lmdzor.calculate_omsk(topo)
    cmsk = lmdzor.calculate_cmsk(omsk)
    assert numpy.all(cmsk == ds_runoff.cmsk.values)