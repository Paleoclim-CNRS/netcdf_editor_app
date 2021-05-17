import numpy
import xarray as xr


def create_ahmcoef(ds_bathy_paleo_orca, bathy_var="Bathymetry"):
    bathy = ds_bathy_paleo_orca[bathy_var].values
    # We stack  bathy horizontally so we don't worry about periodicity
    bathy = numpy.hstack([bathy[:, -6:], bathy, bathy[:, :6]])

    out = numpy.ones(bathy.shape)
    out = out * 100
    out[55:93] = numpy.where(bathy[55:93] > 0, 0, out[55:93])

    # Calculate diffs with periodicity
    def diffs(arr):
        return numpy.vstack([arr[:, -1] - arr[:, 0], numpy.diff(arr).T]).T

    tmp = numpy.ones(out.shape)
    arr = diffs((bathy > 0).astype(int)) * -1
    for i in range(6):
        tmp[arr == -1] = i + 2
        arr = diffs((arr == -1).astype(int))

    out[tmp == 2] = 100
    out[tmp == 3] = 96
    out[tmp == 4] = 69
    out[tmp == 5] = 31
    out[tmp == 6] = 4
    out[tmp == 7] = 0
    out[bathy == 0] = 100
    out[:55] = 100
    out[93:] = 100
    out = out[:, 6:-6]

    dss = xr.Dataset(
        coords={
            "nav_lat": ds_bathy_paleo_orca.nav_lat,
            "nav_lon": ds_bathy_paleo_orca.nav_lon,
        }
    )
    dss["icof"] = ("y", "x"), out
    dss.icof.attrs = {"valid_min": 0.0, "valid_max": 100.0}

    return dss
