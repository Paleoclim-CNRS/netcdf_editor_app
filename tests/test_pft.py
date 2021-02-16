from netcdf_editor_app.utils import pft

import xarray as xr
import numpy

ds_input = xr.open_dataset("./tests/data/input.nc")

def test_generate_pft_netcdf():
    lat_values = numpy.arange(89.75, -90, -0.5)
    lon_vals = numpy.arange(-179.75, 180, 0.5)
    ds = ds_input.interp(
        {
            'lat': lat_values,
            'lon': lon_vals
        }
    )
    topo = ds.topo.values

    latitudes = [15, 35, 50, 80, 90]
    pft_values = numpy.array([
       [ 0, 75, 25,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
       [ 0, 15, 55,  0,  0,  0,  0,  0,  0, 30,  0,  0,  0],
       [ 0,  0,  0,  0, 70, 30,  0,  0,  0,  0,  0,  0,  0],
       [ 0,  0,  0,  0, 40, 30,  0, 30,  0,  0,  0,  0,  0],
       [ 0,  0,  0,  0,  0, 30,  0,  0, 40, 30,  0,  0,  0]])


    ds = pft.generate_pft_netcdf(topo, latitudes, pft_values)
    results = ds.maxvegetfrac.values
    assert results.shape == (1, 13, 360, 720)
    # Make sure ocean values are nan
    x, y = numpy.where(topo < 0)
    assert numpy.all(numpy.isnan(results[0, :, x, y]))

    empty_arrs = [0, 3, 6, 10, 11, 12]
    for index in empty_arrs:
        empty_arr = results[0, index]
        empty_arr = empty_arr[~numpy.isnan(empty_arr)]
        assert numpy.unique(empty_arr) == 0

    # Make sure values are between 0 and 1
    assert numpy.nanmax(results) <= 1
    assert numpy.nanmin(results) >= 0