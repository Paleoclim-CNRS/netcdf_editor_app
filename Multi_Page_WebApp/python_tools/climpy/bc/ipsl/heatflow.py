from datetime import datetime
import numpy
import xarray as xr
from scipy.ndimage import median_filter


def create_heatflow(ds_bathy_paleo_orca, bathy_var="Bathymetry"):
    sfage0 = ((ds_bathy_paleo_orca[bathy_var] - 2600.0) / 365) ** 2
    sfage20 = (-1 / 0.0278) * numpy.log(
        (ds_bathy_paleo_orca.Bathymetry - 5651.0) / (-2473)
    )

    masksf0 = xr.where(sfage0 > 20, 0, 1)
    masksf20 = xr.where(sfage20 >= 20, 1, 0)

    sfage = sfage0 * masksf0 + sfage20 * masksf20

    maskdepth = xr.where(ds_bathy_paleo_orca.Bathymetry > 2500, 1, 0)

    heat0 = 490 / (sfage ** 0.5)
    heat0m = heat0 * maskdepth

    heat0oc = xr.where(heat0m > 400, 400, heat0m)
    heat0ma = xr.where(
        (ds_bathy_paleo_orca.Bathymetry > 0) & (ds_bathy_paleo_orca.Bathymetry < 2500),
        48,
        0,
    )
    heat0ocm = xr.where(heat0oc > 0, heat0oc, 0)
    heat0mam = xr.where(heat0ma > 0, heat0ma, 0)
    heat = heat0ocm + heat0mam

    def xgradient(x):
        return numpy.abs(numpy.gradient(numpy.gradient(x, axis=0), axis=1))

    gradxbathy = xr.apply_ufunc(xgradient, ds_bathy_paleo_orca.Bathymetry)

    maskbathy = xr.where(ds_bathy_paleo_orca.Bathymetry > 0, 1, 0)
    maskgrad = xr.where((gradxbathy > 500) & (heat > 100), 0, 1)

    void_heat = xr.where(heat * maskgrad > 0, heat * maskgrad, numpy.NaN)

    def fill_xy(arr, mask, nb_passes):
        arr = arr.copy(deep=True)
        vals = arr.values
        for _ in range(nb_passes):
            convolution = median_filter(vals, (3, 3))
            vals = numpy.where(mask == 1, convolution, vals)
        arr.values = vals
        return arr

    mask_void = xr.where(
        void_heat.isnull() & (ds_bathy_paleo_orca.Bathymetry > 0), 1, 0
    )
    fill_void_heat = fill_xy(void_heat, mask_void, 2)
    heatflow = fill_void_heat * maskbathy
    ds_out = xr.Dataset()
    ds_out["heatflow"] = heatflow

    # Add coordinates
    ds_out["nav_lon"] = ds_bathy_paleo_orca["nav_lon"]
    ds_out["nav_lat"] = ds_bathy_paleo_orca["nav_lat"]

    # Add time_counter
    ds_model = xr.open_dataset("geothermal_heatingLMParaT.nc", decode_times=False)
    # This adds the time_steps variable and time_counter dim
    ds_out["time_counter"] = ds_model["time_counter"]
    ds_out["time_steps"] = ds_model["time_steps"]
    # Set time_counter on heatflow
    ds_out["heatflow"] = ds_out["heatflow"].expand_dims("time_counter")

    # Set attributes
    ds_out.attrs = {
        'Conventions': 'GDT 1.2',
        'TimeStamp': datetime.now(),
    }
    ds_out["heatflow"].attrs = {
        "units": "mW/m2",
        "title": "heatflow",
        "long_name": "heatflow",
        "Minvalue=": 0.0,
        "Maxvalue=": 200.0,
    }

    return ds_out
