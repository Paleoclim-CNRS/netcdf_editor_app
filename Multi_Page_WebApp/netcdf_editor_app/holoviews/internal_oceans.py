from value_changer import ValueChanger
import panel as pn
import xarray as xr
import numpy

from scipy.ndimage import measurements

import holoviews as hv

colormaps = hv.plotting.list_cmaps()


class InternalOceans(ValueChanger):
    file_type = "raw"
    step = 'internal_oceans'
    elevation_positif = True

    def _calculate_internal_oceans(self):
        # Calculate a binary array of above and below see level
        # from scipy doc:  Any non-zero values in `input` are
        # counted as features and zero values are considered the background.
        # This is why we choose ocean = True
        if self.elevation_positif:
            ocean = self.ds[self.attribute.value] <= 0
        else:
            ocean = self.ds[self.attribute.value] > 0

        # Use scipy to calculate internal oceans
        labeled_array, num_features = measurements.label(ocean)

        # We need to fix the longitude periodicity
        # Algo
        # 1. Find where different
        # 2. Find all connected values
        # 3. Replace par min value
        ids = labeled_array[:, 0] != labeled_array[:, -1]

        # Loop Through and add connected labels
        connections = []
        for val in numpy.where(ids)[0]:
            # Get the values on the edges
            left, right = labeled_array[val, [0, -1]]
            # If we are next to land then we don't care so carry on
            if left == 0 or right == 0:
                continue
            # Variable to test if the values have been encountered before
            added = False
            # loop over each connection
            for i in range(len(connections)):
                # See if one of the values has already been seen
                if left in connections[i] or right in connections[i]:
                    # If it has been seen add both values
                    # TODO probably better way than adding then removing but the arrays shouldnt be too big
                    connections[i].extend([left, right])
                    # Remove the excess values
                    connections[i] = list(set(connections[i]))
                    # Tell the algo not to create a new connections set
                    added = True
                    break
            # If neither of the values have been encountered before then add them
            if not added:
                connections.append([left, right])

        # Replace the values in the array with the minimum seen connected value
        for conn in connections:
            # Get the smallest seen value
            replacement_val = min(conn)
            # For each value replace all the values in the array with the minimum
            # TODO this is will also set the smallest value -> unnecessary
            for val in conn:
                labeled_array[labeled_array == val] = replacement_val

        # Replace continents with numpy.NaN
        # Originally they are ints or floats and numpy.NaN can't be set
        labeled_array = labeled_array.astype(object)
        # continents have a value of 0
        labeled_array[labeled_array == 0] = numpy.NaN
        return labeled_array

    @pn.depends("ds", "attribute.value")
    def load_internal_oceans(self):
        internal_oceans = self._calculate_internal_oceans()
        number_oceans = numpy.nanmax(internal_oceans)

        # Lets counts the number of times each ocean appears this can then be used to
        # Filter out and find the bigger oceans
        nbs, counts = numpy.unique(
            internal_oceans[~numpy.isnan(internal_oceans.astype(float))],
            return_counts=True,
        )

        # Replace the biggest body of water with -1 this will show it as the default body of water
        internal_oceans[internal_oceans == nbs[numpy.argmax(counts)]] = -1

        # Make sure the array shapes line up
        coordinates_shapes = tuple(self.ds.coords.dims.values())
        if internal_oceans.shape == coordinates_shapes:
            internal_oceans = xr.DataArray(internal_oceans, self.ds.coords)
        elif internal_oceans.T.shape == coordinates_shapes:
            internal_oceans = xr.DataArray(internal_oceans.T, self.ds.coords)
        else:
            raise ValueError("Unknown array size of passage problem")

        internal_oceans_image = hv.Image(
            internal_oceans,
            [*self._get_ordered_coordinate_dimension_names()],
            group="Internal_Oceans",
            label=f"Number Internal Oceans: {number_oceans - 1}",
        )
        return internal_oceans_image

    def _get_graphs(self):
        default_grpahs = super()._get_graphs()
        internal_oceans = hv.DynamicMap(self.load_internal_oceans).opts(
            hv.opts.Image(
                "Internal_Oceans",
                clipping_colors={"NaN": "#dedede", "max": "red", "min": "#ffffff"},
                clim=(1.2, 1.5),
                colorbar=False,
                tools=[],
            )
        )
        return default_grpahs + internal_oceans


if "bokeh_app" in __name__:
    int_oc = InternalOceans()
    int_oc.plot().servable("NetCDF Editor")
