from internal_oceans import InternalOceans
import panel as pn
import xarray as xr
import numpy

from scipy.signal import convolve2d

import holoviews as hv

colormaps = hv.plotting.list_cmaps()


class PassageProblems(InternalOceans):
    file_type = "bathy"
    step = "passage_problems"
    elevation_positif = False

    def __init__(self, **params):
        super().__init__(*params)
        self.attribute.value = "topo"
        self.colormap.value = "Batlow"

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
        values = (self.ds[self.attribute.value].values > 0).astype(int)
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

        # potential_points = numpy.zeros(self.ds.trip.values.shape)

        return potential_points

    @pn.depends("ds", "attribute.value")
    def load_passage_problems(self):
        passage_problems = self._calculate_passage_problems()
        number_passage_problems = numpy.sum(passage_problems[passage_problems == 2])

        # Make sure the array shapes line up
        coordinates_shapes = tuple(self.ds.coords.dims.values())
        if passage_problems.shape == coordinates_shapes:
            passage_problems = xr.DataArray(passage_problems, self.ds.coords)
        elif passage_problems.T.shape == coordinates_shapes:
            passage_problems = xr.DataArray(passage_problems.T, self.ds.coords)
        else:
            raise ValueError("Unknown array size of passage problem")

        passage_problems_image = hv.Image(
            passage_problems,
            [*self._get_ordered_coordinate_dimension_names()],
            group="Passage_problems",
            label=f"Number Diffusive Passage cells: {number_passage_problems}",
        )
        return passage_problems_image

    def _update_clims(self):
        # Color clipping occurs at min val NOT included
        # When first setting up viewer for passage problems we do not want to see the bathy = 0 -> the land
        min_value = float(self.ds[self.attribute.value].min()) + 0.5
        max_value = float(self.ds[self.attribute.value].max())
        # Update the limits of the range slider witht the new values
        self.colormap_range_slider.start = min_value
        self.colormap_range_slider.end = max_value
        # Don't necessarily update the min / max values of the colormap
        if self._auto_update_cmap_min:
            self.colormap_min.value = min_value
        if self._auto_update_cmap_max:
            self.colormap_max.value = max_value

    def _get_graphs(self):
        default_grpahs = super()._get_graphs()
        passage_problems = hv.DynamicMap(self.load_passage_problems).opts(
            hv.opts.Image(
                "Passage_problems",
                clipping_colors={"NaN": "#dedede", "max": "red", "min": "#ffffff"},
                clim=(1.2, 1.5),
                colorbar=False,
                tools=[],
            )
        )
        default_grpahs.opts(
            hv.opts.Image(
                "Map",
                clipping_colors={"min": "#dedede", "max": "#ffffff"},
            )
        )
        return default_grpahs + passage_problems


if "bokeh_app" in __name__:
    pp = PassageProblems()
    pp.plot().servable("NetCDF Editor")
