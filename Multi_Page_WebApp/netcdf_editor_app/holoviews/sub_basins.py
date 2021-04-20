from value_changer import ValueChanger
import panel as pn
import xarray as xr
import numpy

import holoviews as hv

from netcdf_editor_app.db import load_file, save_revision, save_step
from netcdf_editor_app.message_broker import send_preprocessing_message

colormaps = hv.plotting.list_cmaps()


class SubBasins(ValueChanger):
    file_type = "bathy"
    step = "subbasins"

    def _default_ocean_values(self, bathy):
        # Atlantic = 1
        # Pacific = 2
        # Indian = 3
        arr = numpy.zeros(bathy.shape)
        arr[:, 0:93] = 2
        arr[:, 93:145] = 1
        arr[:, 145:183] = 3
        arr[bathy <= 0] = 0
        return arr

    def _load_ds(self, _id):
        self.loaded = False
        with self.app.app_context():
            ds = load_file(_id, self.file_type)
        assert ds.Bathymetry.shape == (149, 182)
        ds.Bathymetry.values = self._default_ocean_values(ds.Bathymetry.values)
        ds = ds.rename({"Bathymetry": "Oceans"})
        self.curvilinear_coordinates = None

        number_coordinates_in_system = len(list(ds.coords.variables.values())[0].dims)
        # Standard Grid
        if number_coordinates_in_system == 1:
            pass
        # Curvilinear coordinates
        elif number_coordinates_in_system == 2:
            dims = list(ds[list(ds.coords)[0]].dims)
            # Store the true coordinates for export
            self.curvilinear_coordinates = list(ds.coords)
            # Add the dimension into the coordinates this results in an ij indexing
            ds.coords[dims[0]] = ds[dims[0]]
            ds.coords[dims[1]] = ds[dims[1]]
            # Remove the curvilinear coordinates from the original coordinates
            ds = ds.reset_coords()
        else:
            raise ValueError("Unknown number of Coordinates")
        self.ds = ds
        attributes = list(ds.keys())
        self.attribute.options = attributes
        self.attribute.value = attributes[0]
        self._original_ds = ds.copy(deep=True)
        self.loaded = True
        return True

    def _options_pane_setup(self):
        self.options_pane.clear()
        self.options_pane.extend(
            [
                pn.pane.Markdown("""### Colormaps"""),
                pn.Column(
                    self.colormap,
                    pn.Column(
                        pn.Row(
                            self.colormap_min, pn.layout.HSpacer(), self.colormap_max
                        ),
                        self.colormap_range_slider,
                    ),
                    self.colormap_delta,
                ),
                pn.pane.Markdown("""### Change Values"""),
                pn.Column(
                    self.spinner,
                    self.apply,
                    pn.Row(self.undo_button, self.redo_button),
                ),
            ]
        )

    def _set_values(self, value, calculation_type, selection_expr):
        hvds = hv.Dataset(
            self.ds.to_dataframe(
                dim_order=[*list(self.ds[self.attribute.value].dims)]
            ).reset_index()
        )
        land_indexs = hvds[self.attribute.value] == 0
        hvds.data[self.attribute.value].loc[
            hvds.select(selection_expr).data.index
        ] = value
        hvds.data[self.attribute.value].loc[
            hvds.select(selection_expr).data.index
        ] = value
        hvds.data[self.attribute.value].loc[land_indexs] = 0
        self.ds[self.attribute.value] = tuple(
            (
                list(self.ds[self.attribute.value].dims),
                hvds.data[self.attribute.value].values.reshape(
                    *self.ds[self.attribute.value].shape
                ),
            )
        )
        ds = self.ds.copy(deep=True)
        self.ds = ds

    def _get_graphs(self):
        default_graphs = super()._get_graphs()
        self.colormap.value = "viridis"
        self.colormap_delta.value = 0.75
        # Only allow values 1 to 3
        self.spinner.value = 1
        self.spinner.start = 1
        self.spinner.end = 3
        return default_graphs

    def save(self, event):
        atlmsk = numpy.where(self.ds.Oceans.values == 1, 1, 0)
        pacmsk = numpy.where(self.ds.Oceans.values == 2, 1, 0)
        indmsk = numpy.where(self.ds.Oceans.values == 3, 1, 0)
        ds = xr.Dataset(
            coords={},
            data_vars={
                "navlat": (["y", "x"], self.ds.nav_lat.values),
                "navlon": (["y", "x"], self.ds.nav_lon.values),
                "atlmsk": (["y", "x"], atlmsk),
                "pacmsk": (["y", "x"], pacmsk),
                "indmsk": (["y", "x"], indmsk),
            },
        )

        ds["navlon"].attrs = {"units": "degrees_east"}
        ds["navlat"].attrs = {"units": "degrees_north"}
        ds["atlmsk"].attrs = {}
        ds["pacmsk"].attrs = {}
        ds["indmsk"].attrs = {}

        with self.app.app_context():
            save_revision(self.data_file_id, ds, "sub_basins")
            if self.step is not None:
                save_step(
                    self.data_file_id,
                    step=self.step,
                    parameters={"id": self.data_file_id},
                    up_to_date=True,
                )
                send_preprocessing_message(
                    self.step + ".done", message={"id": self.data_file_id}
                )


if "bokeh_app" in __name__:
    sub_basins = SubBasins()
    sub_basins.plot().servable("NetCDF Editor")
