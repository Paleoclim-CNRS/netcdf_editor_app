from netcdf_editor_app import create_app
from netcdf_editor_app.db import load_file, save_revision
from bokeh.models import FixedTicker
import panel as pn
from holoviews import opts
from holoviews.selection import link_selections
import xarray as xr
import param
import numpy

from scipy.ndimage import measurements
from scipy.signal import convolve2d

import hvplot.xarray
import holoviews as hv

import palettable


CUSTOM_COLORMAPS = {
    'Oleron': palettable.scientific.sequential.Oleron_18.mpl_colormap,
    'Vik': palettable.scientific.diverging.Vik_18.mpl_colormap,
    'Cork': palettable.scientific.diverging.Cork_18.mpl_colormap,
    'LaPaz': palettable.scientific.sequential.LaPaz_18.mpl_colormap,
    'Batlow': palettable.scientific.sequential.Batlow_18.mpl_colormap
}

colormaps = hv.plotting.list_cmaps()
colormaps.extend(CUSTOM_COLORMAPS.keys())
colormaps = sorted(colormaps, key=lambda L: (L.lower(), L))

opts.defaults(
    opts.Image(
        # Values taken from holoviews.Store.custom_options for a xarray.Dataset.hvplot()
        colorbar=True,
        logx=False,
        logy=False,
        responsive=True,
        aspect=2,
        shared_axes=True,
        show_grid=False,
        show_legend=True,
        tools=["hover", "lasso_select", "box_select"],  # Default = hover
    )
)

pn.config.sizing_mode = "stretch_width"


class ValueChanger(param.Parameterized):

    # Holoviews.DataSet => Data
    ds = param.Parameter()
    # Used to store when inital data is loaded
    loaded = param.Parameter()

    def __init__(self, **params):
        # How we are going to modify the values
        # Absolute => Set to that value
        # Relatif => Base value + new value
        # Percentage => Base value + percentage
        self.calculation_type = pn.widgets.RadioButtonGroup(
            options=["Absolute", "Relatif", "Percentage"], align="end"
        )
        # Replacement value
        self.spinner = pn.widgets.IntInput(
            name="Replacement Value", value=0, align="start"
        )

        # Buttons
        self.apply = pn.widgets.Button(
            name="\u2713 Apply", align="end", button_type="primary"
        )
        self.undo_button = pn.widgets.Button(
            name="\u21B6 Undo", align="end", button_type="warning"
        )
        self.redo_button = pn.widgets.Button(
            name="Redo \u21B7", align="end", button_type="warning"
        )
        self.save_button = pn.widgets.Button(
            name="Save", align="end", button_type="success"
        )
        # Mask
        self.mask = pn.widgets.Checkbox(name="Mask", max_width=100)
        self.mask_value = pn.widgets.IntInput(name="Mask Value", value=0)

        # Store the variable we want to look at and modify
        self.attribute = pn.widgets.Select(name="Variable", max_width=200, align="end")
        # Choose colormap
        self.colormap = pn.widgets.Select(
            name="Colormap",
            options=colormaps,
            value='Oleron',
            max_width=200,
            align="start",
        )
        self.colormap_min = pn.widgets.IntInput(name="Min Value", width=100)
        self.colormap_max = pn.widgets.IntInput(
            name="Max Value", width=100, align="end"
        )
        self.colormap_range_slider = pn.widgets.RangeSlider(width=400, show_value=False)
        self.colormap_delta = pn.widgets.FloatInput(
            name="Delta between values", value=0, align="end"
        )
        # Link the viewing of multiple graphs together
        self.selection = link_selections.instance(unselected_alpha=0.4)

        # Parts of the display
        self.file_pane = pn.Column()
        self.graph_pane = pn.Column()
        self.options_pane = pn.Column()

        self.param.ds.default = xr.Dataset()
        self.param.loaded.default = False
        super().__init__(**params)

        self.apply.on_click(self._apply_values)
        self.undo_button.on_click(self.undo)
        self.redo_button.on_click(self.redo)
        self.save_button.on_click(self.save)
        self.save_button.js_on_click(
            args={
                "target": pn.state.curdoc.session_context.request.arguments["redirect"][
                    0
                ].decode()
            },
            code="window.location.href = target",
        )
        self._auto_update_cmap_min = True
        self._auto_update_cmap_max = True

        self.curvilinear_coordinates = None
        self._undo_list = []
        self._redo_list = []

        self.colormap_min.param.watch(self._colormap_callback, "value")
        self.colormap_max.param.watch(self._colormap_callback, "value")
        self.colormap_range_slider.param.watch(self._colormap_callback, "value")
        self.app = create_app()
        try:
            self.file_type = pn.state.curdoc.session_context.request.arguments["file_type"][0].decode()
        except KeyError:
            pass
        self.data_file_id = int(
            pn.state.curdoc.session_context.request.arguments["id"][0]
        )
        self._load_ds(self.data_file_id)
        self._options_pane_setup()

    def _load_ds(self, _id):
        self.loaded = False
        with self.app.app_context():
            ds = load_file(_id, self.file_type)
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

    def _colormap_callback(self, *events):
        event = events[0]
        if event.obj == self.colormap_min:
            # The colormap value has been changed.
            # We need to potentially update the colormap range slider min value
            # We need to potentially update the colormap_range_slider values
            new_value = int(event.new)
            if self.colormap_range_slider.start > new_value:
                self.colormap_range_slider.start = new_value

            vals = list(self.colormap_range_slider.value)
            if vals[0] != new_value:
                vals[0] = new_value
                self.colormap_range_slider.value = tuple(vals)

            # Have we manually changed the minimum value?
            min_value = int(self.ds[self.attribute.value].min())

            if new_value == min_value:
                # The new value is the minimum value so we want the min values to
                # update automatically
                self._auto_update_cmap_min = True
            else:
                self._auto_update_cmap_min = False

        elif event.obj == self.colormap_max:
            # The colormap value has been changed.
            # We need to potentially update the colormap range slider max value
            # We need to potentially update the colormap_range_slider values
            new_value = int(event.new)
            if self.colormap_range_slider.end < new_value:
                self.colormap_range_slider.end = new_value

            vals = list(self.colormap_range_slider.value)
            if vals[1] != new_value:
                vals[1] = new_value
                self.colormap_range_slider.value = tuple(vals)

            # Have we manually changed the max value?
            max_value = int(self.ds[self.attribute.value].max())
            if new_value == max_value:
                self._auto_update_cmap_max = True
            else:
                self._auto_update_cmap_max = False

        elif event.obj == self.colormap_range_slider:
            # Lets see whcih values have changed
            new_vals = event.new
            old_vals = event.old

            # Minimum value has changed
            if new_vals[0] != old_vals[0]:
                self.colormap_min.value = int(new_vals[0])
            # Maximum value has changed
            if new_vals[1] != old_vals[1]:
                self.colormap_max.value = int(new_vals[1])

    def _set_values(self, value, calculation_type, selection_expr):
        hvds = hv.Dataset(
            self.ds.to_dataframe(
                dim_order=[*list(self.ds[self.attribute.value].dims)]
            ).reset_index()
        )
        if calculation_type == "Absolute":
            hvds.data[self.attribute.value].loc[
                hvds.select(selection_expr).data.index
            ] = value
        elif calculation_type == "Relatif":
            hvds.data[self.attribute.value].loc[
                hvds.select(selection_expr).data.index
            ] += value
        elif calculation_type == "Percentage":
            hvds.data[self.attribute.value].loc[
                hvds.select(selection_expr).data.index
            ] *= (100 + value) / 100.0
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

    def undo(self, event):
        # Nothing in the undo list
        if not len(self._undo_list):
            return

        # Get the last action in the undo list
        undo_action = self._undo_list.pop()
        self._redo_list.append(undo_action.copy())

        # If it is 'Absolute' Change we don't stock the
        # initial values so we have to run all the steps up to this one to
        # undo this change
        if undo_action["calculation_type"] in ["Absolute", "Depression_filling"]:
            # We reset the dataset to it's initial value
            self.ds = self._original_ds.copy(deep=True)
            # We apply each step one by one
            for action in self._undo_list:
                self._apply_action(action)

        elif undo_action["calculation_type"] == "Relatif":
            # Apply the opposite transformation
            undo_action["value"] *= -1
            self._apply_action(undo_action)
        elif undo_action["calculation_type"] == "Percentage":
            # Apply the opposite transformation
            undo_action["value"] = ((100 * 100) / (100 + undo_action["value"])) - 100
            self._apply_action(undo_action)
        else:
            raise ValueError(
                "Can not undo action, unknown calculation type {}".format(
                    undo_action["calculation_type"]
                )
            )

    def redo(self, event):
        # Nothing in the redo list
        if not len(self._redo_list):
            return

        # Get the last action in the redo list
        redo_action = self._redo_list.pop()

        self._apply_action(redo_action)
        # Add the action to the list of undo actions
        self._undo_list.append(redo_action)

    def save(self, event):
        with self.app.app_context():
            save_revision(self.data_file_id, self.ds, self.file_type)

    def _apply_action(self, action):
        if action["calculation_type"] in ["Absolute", "Percentage", "Relatif"]:
            self._set_values(
                value=action["value"],
                calculation_type=action["calculation_type"],
                selection_expr=action["selection_expr"],
            )
        else:
            raise ValueError(
                "Cannot apply step {}, unknown calculation_type {}".format(
                    action, action["calculation_step"]
                )
            )

    def _apply_values(self, event):
        if self.selection.selection_expr is None:
            return
        action = {
            "selection_expr": self.selection.selection_expr,
            "calculation_type": self.calculation_type.value,
            "value": self.spinner.value,
        }
        # Apply the action
        self._apply_action(action)

        # Add the action to the list of undo actions
        self._undo_list.append(action)

        self.selection.selection_expr = None

    def _get_ordered_coordinate_dimension_names(self):
        dimension_names = list(self.ds.coords)
        if "lat" in dimension_names[0].lower() and "lon" in dimension_names[1].lower():
            dimension_names = dimension_names[::-1]
        elif "x" == dimension_names[1].lower() or "y" == dimension_names[0].lower():
            dimension_names = dimension_names[::-1]
        return dimension_names

    def _calculate_internal_oceans(self):
        # Calculate a binary array of above and below see level
        # from scipy doc:  Any non-zero values in `input` are
        # counted as features and zero values are considered the background.
        # This is why we choose ocean = True
        ocean = self.ds[self.attribute.value] <= 0

        # Use scipy to calculate internal oceans
        labeled_array, num_features = measurements.label(ocean)

        # Replace continents with numpy.NaN
        # Originally they are ints or floats and numpy.NaN can't be set
        labeled_array = labeled_array.astype(object)
        # continents have a value of 0
        labeled_array[labeled_array == 0] = numpy.NaN

        return labeled_array

    def _calculate_passage_problems(self):
        # Define template we are looking for passages
        # Where only diffusion occurs this means we are looking
        # for ocean passages one in width/height
        # 1 => Ocean
        # -1 => Land
        # 0 = Indifferent
        template = numpy.array([[0, 1, 0], [-1, 1, -1], [0, 1, 0]])

        # Theoretical max value when the template is found
        # Note that 0s are considered wildcards so they are not taken into
        # Account
        # TODO this only works on data arrays where the absolute values are 1
        perfect_match = numpy.sum(numpy.abs(template))

        # we recode the values of land to -1 as
        # we did in the template
        values = (self.ds[self.attribute.value].values <= 0).astype(int)
        values[values == 0] = -1

        # Create an empty array where we are going to stock the values
        # TODO This could potentially by a binary array??
        potential_points = values
        #         potential_points[:] = numpy.nan

        # Mark points where there is only diffusion in longitude direction
        convolvedh = convolve2d(values, template, "same")
        potential_points[convolvedh == perfect_match] = 2

        # Mark points where there is only diffusion in latitude direction
        convolvedv = convolve2d(values, template.T, "same")
        potential_points[convolvedv == perfect_match] = 2

        potential_points = potential_points.astype(object)
        potential_points[potential_points == -1] = numpy.NaN

        return potential_points

    def _options_pane_setup(self):
        self.options_pane.clear()
        self.options_pane.extend(
            [
                pn.pane.Markdown("""### Variable"""),
                pn.Column(self.attribute),
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
                pn.pane.Markdown("""### Mask"""),
                pn.Row(self.mask, self.mask_value),
                pn.pane.Markdown("""### Change Values"""),
                pn.Column(
                    self.calculation_type,
                    self.spinner,
                    self.apply,
                    pn.Row(self.undo_button, self.redo_button),
                ),
            ]
        )

    def get_grid_style(self):
        # Calculate Ticks
        ydim, xdim = self.ds[self.attribute.value].dims
        xvals = self.ds[xdim].values
        yvals = self.ds[ydim].values
        x_ticks = (xvals[1:] + xvals[:-1]) / 2
        y_ticks = (yvals[1:] + yvals[:-1]) / 2
        # Setup a grid style
        grid_style = {
            "grid_line_color": "black",
            "grid_line_width": 1,
            "xgrid_ticker": x_ticks,
            "ygrid_ticker": y_ticks,
        }
        return grid_style

    def _cmap(self):
        if self.colormap.value in CUSTOM_COLORMAPS.keys():
            return CUSTOM_COLORMAPS[self.colormap.value]
        return self.colormap.value

    def _update_clims(self):
        min_value = int(self.ds[self.attribute.value].min())
        max_value = int(self.ds[self.attribute.value].max())
        # Update the limits of the range slider witht the new values
        self.colormap_range_slider.start = min_value
        self.colormap_range_slider.end = max_value
        # Don't necessarily update the min / max values of the colormap
        if self._auto_update_cmap_min:
            self.colormap_min.value = min_value
        if self._auto_update_cmap_max:
            self.colormap_max.value = max_value

    def _clims(self):
        if self.mask.value:
            return self.mask_value.value, self.mask_value.value
        else:
            return self.colormap_min.value, self.colormap_max.value

    def _color_levels(self):
        if self.colormap_delta.value <= 0:
            return None
        return (
            list(
                numpy.arange(
                    self.colormap_min.value,
                    self.colormap_max.value,
                    self.colormap_delta.value,
                )
            )
            + [self.colormap_max.value]
        )

    def _colorbar_opts(self):
        if self.colormap_delta.value <= 0:
            return {}
        ticks = self._color_levels()
        if len(ticks) > 8:
            ticks = ticks[:: len(ticks) // 8] + [ticks[-1]]
        # Add 0 to the ticks
        if self.colormap_min.value * self.colormap_max.value < 0:  # Either side of 0
            ticks = numpy.insert(ticks, numpy.searchsorted(ticks, 0), 0)
        return {"ticker": FixedTicker(ticks=ticks)}

    @pn.depends(
        "colormap.value",
        "colormap_min.value",
        "colormap_max.value",
        "mask.value",
        "mask_value.value",
        "colormap_delta.value",
    )
    def _opts(self, element):
        return element.opts(
            cmap=self._cmap(),
            clim=self._clims(),
            color_levels=self._color_levels(),
            colorbar_opts=self._colorbar_opts(),
        )

    @pn.depends("ds", "attribute.value")
    def load_attribute_map(self):
        self._update_clims()
        return hv.Image(
            self.ds[self.attribute.value],
            [*self._get_ordered_coordinate_dimension_names()],
        )

    def _get_graphs(self):
        return hv.DynamicMap(self.load_attribute_map).apply(self._opts)

    @pn.depends("loaded", watch=True)
    def get_plots(self):
        if not self.loaded:
            return

        graphs = self._get_graphs()

        layout = self.selection(graphs + self.ds[self.attribute.value].hvplot.hist())

        layout.opts(
            hv.opts.Histogram(tools=["hover"]),
            hv.opts.Image(
                tools=["hover", "box_select", "lasso_select"],
                show_grid=True,
                gridstyle=self.get_grid_style(),
                alpha=0.75,
            ),
        ).cols(1)

        self.graph_pane.clear()
        self.graph_pane.append(layout)
        self._auto_update_cmap_min = True
        self._auto_update_cmap_max = True

    def plot(self):
        template = pn.template.MaterialTemplate(
            title="NetCDF Editor App",
            logo="https://raw.githubusercontent.com/CEREGE-CL/CEREGE-CL.github.io/main/logo.png",
            favicon="https://raw.githubusercontent.com/CEREGE-CL/CEREGE-CL.github.io/main/logo.png",
            header_background="#42a5f5",
        )
        template.sidebar.append(self.file_pane)
        template.sidebar.append(self.options_pane)
        template.main.append(self.graph_pane)
        template.main.append(self.save_button)
        return template


if "bokeh_app" in __name__:
    vc = ValueChanger()
    vc.plot().servable("NetCDF Editor")
