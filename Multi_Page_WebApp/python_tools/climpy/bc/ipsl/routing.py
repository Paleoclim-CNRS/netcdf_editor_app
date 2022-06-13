from datetime import datetime
import numpy
import os

from scipy.signal import convolve2d
from scipy.ndimage import measurements
from skimage.morphology import reconstruction

import xarray as xr

next_downstream_cell = numpy.array(
    [
        # 0  1  2  3   4   5   6   7   8
        [0, 1, 1, 0, -1, -1, -1, 0, 1],
        [0, 0, 1, 1, 1, 0, -1, -1, -1],
    ]
).T

# Set the trip values, this is in case we don't use a standard order in arrays
trip_values = numpy.array([1, 2, 3, 4, 5, 6, 7, 8])


def get_adjacent_coords(i, j, arr):
    imax, jmax = arr.shape
    if j == jmax - 1:
        jp1 = 0
    else:
        jp1 = j + 1
    if j == 0:
        jm1 = jmax - 1
    else:
        jm1 = j - 1

    ip1 = numpy.min((i + 1, imax - 1))
    im1 = numpy.max((0, i - 1))

    return (
        numpy.array(numpy.meshgrid([im1, i, ip1], [jm1, j, jp1], indexing="ij"))
        .reshape(2, -1)
        .T
    )


def get_adjacent_values(i, j, arr):
    return numpy.array(
        [arr[ii, jj] for ii, jj in get_adjacent_coords(i, j, arr)]
    ).reshape(3, -1)


def get_padded_array(arr, ensure_gradient=False):
    # Add a pad of 1 around the edge and use the wrap method
    padded = numpy.pad(arr, 1, "wrap")
    # Fix the Top and bottom because they are the same value as before and not wrapped (there is no North/South wrapping)
    padded[0] = padded[1]
    padded[-1] = padded[-2]
    if ensure_gradient:
        padded[0] += 1.0
        padded[-1] += 1
    return padded


def _fill_depressions_in_topo(dem):
    # Copy the values to make sure we are not editing the original values
    dem = numpy.copy(dem)
    # Set Antartic values (row 0) with the maximum value of each column for the first two rows
    # This should ensure that water will never run towards to pole
    # We increment the value slightly to ensure it is not flat
    dem[0] = numpy.max(dem[:2], axis=0) + 1.0
    # Copy the map to either side for east - west periodicity
    # Rather than hard coding the periodicity as the array is small and
    # the function runs quickly we fix the topography aligned 3 times next to each other
    dem = numpy.hstack((dem, dem, dem))

    # Create a pad of one all the way around the dem and set to max value
    # Scipy reconstruction works bby eroding awy the values
    # See here:
    # https://scikit-image.org/docs/dev/auto_examples/features_detection/plot_holes_and_peaks.html#sphx-glr-auto-examples-features-detection-plot-holes-and-peaks-py
    padded = numpy.pad(dem, 1, "maximum")
    seed = numpy.copy(padded)
    # Add borders
    seed[1:-1, 1:-1] = dem.max() + 1

    filled = reconstruction(seed, padded, method="erosion")
    # Because we padded the array we need to remove this
    filled = filled[1:-1, 1:-1]
    # Reset the ocean -> dont modify the ocean
    filled[dem <= 0] = dem[dem <= 0]
    # Take the middle array (we stacked 3 of them next to each other)
    filled = filled[:, int(filled.shape[1] / 3) : 2 * int(filled.shape[1] / 3)]
    return filled


def _migrate_front_across_flats(front, arrays, trip_values, iterations, iteration):
    """
    Recursive function to propagate values across a front. We use recursion rather than a while loop
    """
    # Get the adjacent cells that have the same value as the center cell
    # arrays is defined below and is a 8 by i by j where 8 corresponds the the different directions
    # each layer of arrays corresponds to the difference in the trip value
    # in that direction eg layer one is the differences with cells to the north
    flat_coords = numpy.argwhere(arrays[:, front[0], front[1]].T == 0)
    # if there are no more cells with adjacent flat cells then stop the recursion
    if not len(flat_coords):
        return

    # Get the adjacent cell coordinates by using where the flat cell is in comparison to the each cell on the front
    # Downstream cells points to the combination of -1, 0, 1 in each direction that needs adding
    next_cells = (
        next_downstream_cell[trip_values[flat_coords[:, 1]]]
        + front.T[flat_coords[:, 0]]
    )
    # Fix the EW periodicity
    # Replace values that are -1 with the last value
    next_cells[:, 1] = numpy.where(
        next_cells[:, 1] == -1, iterations.shape[1] - 1, next_cells[:, 1]
    )
    # Replace values that go off the end with 0
    next_cells[:, 1] = numpy.where(
        next_cells[:, 1] == iterations.shape[1], 0, next_cells[:, 1]
    )

    # Remove duplicate cells that can be found from multiple routes
    # This takes the first occurence each time -> means that the trip values aren't going to be random
    # They are probably ordered in the following order 6, 5, 4, 7, 0, 3, 8, 1, 2 (trip values order)
    next_cells, unique_ids = numpy.unique(next_cells, axis=0, return_index=True)
    # Add the trip values to the next cells this is for consitent filtering
    next_cells = numpy.hstack(
        (next_cells, trip_values[flat_coords[unique_ids, 1]].reshape(-1, 1))
    )
    # Make sure we don't have any values that go outside the array
    # TODO this may be casuing problems ??? but hopefully isn't being hit
    next_cells = next_cells[
        (next_cells[:, 0] >= 0)
        & (next_cells[:, 0] < iterations.shape[0])
        & (next_cells[:, 1] >= 0)
        & (next_cells[:, 1] < iterations.shape[1])
    ]
    # only take values that haven't already been seen
    # We do this by checking if a value for the cell has been set ->
    # this means it has been seen in another layer of the recursion
    # trip values are 1-> 8 this means that anythin 0 or below hasn't been seen
    next_cells = next_cells[iterations[next_cells[:, 0], next_cells[:, 1]] <= 0]
    # Change the values of the next cells to their correct trip values (calculated where they came from)
    iterations[next_cells[:, 0], next_cells[:, 1]] = iteration
    # Redo this function as many times as necessary
    iteration += 1
    return _migrate_front_across_flats(
        next_cells[:, :2].T, arrays, trip_values, iterations, iteration
    )


def get_differences_with_neighbors(array, ensure_gradient=False):
    # Pad the array so that we don't worry about EW periodicity
    arr = get_padded_array(array, ensure_gradient=ensure_gradient)
    # Calculate columnwise differences -> NS differences
    NS = arr[1:] - arr[:-1]
    # Calculate rowwise differnces -> EW differences
    EW = arr[:, :-1] - arr[:, 1:]
    # Calculate diagonals
    NESW = arr[1:, 1:] - arr[:-1, :-1]
    SENW = arr[1:, :-1] - arr[:-1, 1:]

    # Put it all into a multidimensional array
    # Because we calculated differences in different directions we need to take parts of each array
    # Also because we calculated the difference in one direction the value in the other direction is the symmetric (*-1)
    # Big positive values means the biggest differences -> direction of descent
    arrays = numpy.array(
        [
            NS[1:, 1:-1] * -1,
            NESW[1:, 1:] * -1,
            EW[1:-1, 1:],
            SENW[:-1, 1:],
            NS[:-1, 1:-1],
            NESW[:-1, :-1],
            EW[1:-1, :-1] * -1,
            SENW[1:, :-1] * -1,
        ]
    )
    return arrays


def _add_gradient_to_flats(topo):
    """
    calculate trip values using matrices rather than loops
    """
    arrays = get_differences_with_neighbors(topo, ensure_gradient=True)

    # Fix the values on flats
    # There is at least one cell where the water can flow downwards
    # pas is like in the mountains not quite a col but like a suspended glacial valley
    is_pas = numpy.sum(arrays > 0, axis=0) > 0
    # A pas has a least one value that is greater than 0 -> downwards flow
    # Flats are defined by areas where at least the difference between 2 cells is 0
    is_flat = numpy.sum(arrays == 0, axis=0) > 0

    # Only take land cells
    # TODO maybe we should be passing omsk in to ensure we have the correct values
    is_pas[topo <= 0] = False
    is_flat[topo <= 0] = False

    # Exit points from flat basins are defined as the points that have the same
    # altitude as a neighbor cell (is_flat) and at least one downward cell (is_pas)
    exit_points = numpy.array(
        numpy.where((is_pas == True) & (is_flat == True))  # noqa: E712
    )

    # Store an array of when the cell was seen for the first time
    # This is passed to our recursively function and filled over time
    iterations = numpy.zeros(topo.shape)
    iteration = 1

    # Add the values to iterations this is so the algo knows that these
    # cells ahev already been seen and will expand / migrate from here
    iterations[exit_points[0], exit_points[1]] = iteration
    iteration += 1

    # Migrate the front across the flats recursively
    _migrate_front_across_flats(exit_points, arrays, trip_values, iterations, iteration)
    zeros = numpy.where(iterations == 0)
    ones = numpy.where(iterations == 1)
    #     iterations += numpy.random.rand(*iterations.shape) / 100
    iterations[zeros] = 0
    iterations[ones] = 1
    increment = 1.0 / (100 * iterations.max())
    topo += iterations * increment

    return topo


def _ensure_south_up(topo, latitudes):
    if latitudes[0] < 0:
        assert latitudes[-1] > 0
        return topo
    if latitudes[-1] < 0:
        assert latitudes[0] > 0
        return topo[::-1]
    raise AssertionError("Can not find South, latitude values do not vary across 0")


def fix_topo(topo, latitudes):
    topo = topo.copy()
    # Some files with trailing float problems (1e-14) have been observed this should fix that
    topo = numpy.round(topo, 3)
    topo = _ensure_south_up(topo, latitudes)
    topo = _fill_depressions_in_topo(topo)
    topo = _add_gradient_to_flats(topo)
    return topo


def calculate_orog(topo):
    orog = topo.copy()
    # Set see  values to NaN
    # TODO pass omsk in??
    orog[orog <= 0] = numpy.nan
    return orog


def calculate_omsk(topo):
    # omsk has the value of 0 on continents and 1 in ocean
    # Oceans are defined as values less than or equal to 0
    omsk = numpy.zeros(topo.shape)
    omsk[topo <= 0] = 1
    return omsk


def calculate_cmsk(omsk):
    template = numpy.ones(9).reshape(3, 3)
    # pad the array so we don't have to worry about edge cases

    is_water = get_padded_array(omsk).astype(bool)
    # Calculate the convolution
    conv = convolve2d(is_water, template, "same")
    # Extract the values we want
    # Conv is positive where one of the adjacent cells is ocean
    # We only want land cells
    res = (conv > 0) & (~is_water)
    # We padded omsk so as to not worry about the periodicity
    # Now we have to remove the padding
    res = res[1:-1, 1:-1]
    return res


def calculate_continents(topo, omsk):

    # Calculate the continents measurements returns a list of distinct objects in the image
    continents = measurements.label(~omsk.astype(bool))[0]
    # Take a look at the edges and see where continents are different and should be the same due to cyclicite
    difference_long_180 = numpy.argwhere(
        (continents[:, 0] != continents[:, -1])
        &  # Check first column and last column are the same
        # Make sure the first column isn't ocean -> last column is a coast
        (continents[:, 0] != 0)
        &
        # Make sure the last column isn't ocean -> first column is a coast
        (continents[:, -1] != 0)
    ).flatten()
    same_continents = []
    for d in difference_long_180:
        vals = [continents[d, 0], continents[d, -1]]
        inserted = False
        for i in range(len(same_continents)):
            conts = same_continents[i]
            # At least one of the values already exists in a list
            if len(set(vals) & set(conts)):
                conts.extend(vals)
                same_continents[i] = list(set(conts))
                inserted = True
                break
        # If the continents haven't already been seen then add a new continent to the list
        if not inserted:
            same_continents.append(vals)

    for cont in same_continents:
        # Replace all values in continents with the first value
        for i in range(1, len(cont)):
            continents[continents == cont[i]] = cont[0]
    return continents


def calculate_trip(topo, omsk):
    """
    calculate trip values using matrices rather than loops
    """
    arrays = get_differences_with_neighbors(topo, ensure_gradient=True)

    # Find the index of the biggest value -> direction of descent if multiple values occur it takes the smallest (first seen)
    ind = numpy.argmax(arrays, axis=0)
    # Convert indexs to trip values
    rtm = trip_values[ind].astype(float)
    # only take rtm values on land
    rtm[omsk == True] = numpy.nan  # noqa: E712

    # TODO we probably need to ensure topo always has the south Up
    # The first row (Antartic) is considered to be the pole so say everything goes north from there
    rtm[0] = numpy.where(~numpy.isnan(rtm[0]), 1, numpy.NaN)

    # The last row (Artic) is considered to be the pole so say everything goes South from there
    rtm[-1] = numpy.where(~numpy.isnan(rtm[-1]), 5, numpy.NaN)

    return rtm


def calculate_curvilinear_coordinates():
    # Calculate longigtude and latitude array
    rlon = numpy.arange(-179.5, 180.5)
    rlat = numpy.arange(89.5, -90.5, -1)
    return numpy.meshgrid(rlon, rlat)


def _calculate_dx(rlat):
    rEarth = 6370
    return (2 * numpy.pi * rEarth * numpy.cos(rlat * numpy.pi / 180)) / 360.0


def _calculate_dy():
    rEarth = 6370
    return (numpy.pi * rEarth) / 180


def calculate_area(rlat):
    dy = _calculate_dy()
    dx = _calculate_dx(rlat)
    area = dx * dy
    return area


def get_next_cell(cells, temp_array, trip):
    """
    Recursive function that from an given cell will find all the next down stream cells using trip values
    """
    ii, jj = cells[-1]
    # Subscripts pointing on the next cell (downstream) following trip value
    downstream_direction = next_downstream_cell[int(trip[ii, jj])]
    # Add the indexs to the original values
    ip1 = ii + downstream_direction[0]
    jp1 = jj + downstream_direction[1]
    # Handle periodicity
    if jp1 < 0:
        jp1 = temp_array.shape[1] - 1
    if jp1 >= temp_array.shape[1]:
        jp1 = 0
    next_cell = (ip1, jp1)
    if ip1 >= temp_array.shape[0] or ip1 < 0:
        raise AssertionError(f"{datetime.now()} Shouldn't be here, {ii}, {jj}")
        print(f" {datetime.now()} Shouldn't be here, {ii}, {jj}")
        return cells, "outside", next_cell
    # If the next cell is ocean then we are done
    if numpy.isnan(trip[next_cell]):
        return cells, "ocean", next_cell
    # If the next cell already has a value then we are done
    if temp_array[next_cell] != 0:
        return cells, "junction", next_cell

    # We have already come across this cell lets see if we can find a different output
    if next_cell in cells:
        return cells, "infinite_loop", next_cell

    cells.append(next_cell)
    return get_next_cell(cells, temp_array, trip)


def calculate_basins(topo, trip, area):
    # create dummy array
    basins = numpy.zeros(topo.shape)
    # Find the highest point
    highest_point = numpy.unravel_index(
        numpy.where(basins == 0, topo, -9999).argmax(), topo.shape
    )

    # Setup a basin number
    basin_nb = 1

    # While we still have points that are on land (higher than 0)
    #  calculate the runoff cells from this high point
    while topo[highest_point] > 0:
        # Get cells connected to current cell
        cells, end_reason, next_cell = get_next_cell([highest_point], basins, trip)
        x, y = numpy.array(cells).T
        # If we hit a junction then we need to adjust all the values
        if end_reason == "junction":
            basins[x, y] = basins[next_cell]
        else:
            # Set value to a new basin
            basins[x, y] = basin_nb
            basin_nb += 1
        # Calculate the next highest point
        highest_point = numpy.unravel_index(
            numpy.where(basins == 0, topo, -9999).argmax(), topo.shape
        )

    # Rarrange by the area
    # We use bincount with weights corresponding to the area of each cell
    sizes = numpy.bincount(basins.flatten().astype(int), area.flatten())
    # We sort the cells by size
    sorted_sizes = numpy.argsort(sizes)[::-1]
    # We assign the new value
    for i in range(len(sorted_sizes)):
        basins[basins == sorted_sizes[i]] = -i

    # REplace ocean with zeros
    basins[basins == 0] = numpy.nan
    # Flip to get positive numbers
    basins *= -1
    return basins


def calculate_ocean_distances(trip, rlat):
    dx = _calculate_dx(rlat)
    dy = _calculate_dy()
    # Creat Distance for each cell based on trip values
    trip_distances = numpy.array(
        [
            # Dummy used for indexing this should never be choosen
            numpy.zeros(trip.shape),
            numpy.ones(trip.shape) * dy,
            numpy.sqrt(dx ** 2 + dy ** 2),
            dx,
            numpy.sqrt(dx ** 2 + dy ** 2),
            numpy.ones(trip.shape) * dy,
            numpy.sqrt(dx ** 2 + dy ** 2),
            dx,
            numpy.sqrt(dx ** 2 + dy ** 2),
        ]
    )

    temp_trip = numpy.where(~numpy.isnan(trip), trip.astype(int), 0)
    xx, yy = numpy.ix_(numpy.arange(trip.shape[0]), numpy.arange(trip.shape[1]))
    distances = trip_distances[temp_trip, xx, yy]
    distances[distances == 0] = numpy.nan

    return distances


def calculate_outflow_points(basins, trip):
    # We store the outflow points as they are used later to determine trip values
    # They could be calculated using minimum topo height inside each basin but this has problems
    # When the outflow point is at the same height as other points
    # Considering we are following flow directions we will just store the values
    outflow_points = []
    for basin_nb in range(
        numpy.nanmin(basins).astype(int), numpy.nanmax(basins).astype(int) + 1
    ):
        nboutflow = 0
        outloc = [numpy.nan, numpy.nan]

        # Calculate the nb of outflow points
        x, y = numpy.where(basins == basin_nb)
        for k in range(len(x)):
            i, j = x[k], y[k]
            downstream_direction = next_downstream_cell[int(trip[i, j])]
            # Add the indexs to the original values
            ip1 = i + downstream_direction[0]
            jp1 = j + downstream_direction[1]
            # Handle periodicity
            if jp1 < 0:
                jp1 = basins.shape[1] - 1
            if jp1 >= basins.shape[1]:
                jp1 = 0
            next_cell = (ip1, jp1)
            if ~numpy.isnan(trip[i, j]) & numpy.isnan(trip[next_cell]):
                if outloc != [i, j]:
                    outloc = [i, j]
                    nboutflow += 1
        outflow_points.append(outloc)

        if nboutflow > 1:
            raise AssertionError(
                "error occured to many outflow points for basin {}".format(basin_nb)
            )
            print("error occured to many outflow points")
        elif nboutflow == 0:
            raise AssertionError("No outflow points for basin {}".format(basin_nb))
            print("No outflow points for basin {}".format(basin_nb))

    outflow_points = numpy.array(outflow_points)
    return outflow_points


def calculate_river_lengths(topo, trip, ocean_distances, omsk):
    # create dummy array
    flength = numpy.zeros(trip.shape)
    # Find the highest point
    highest_point = numpy.unravel_index(
        numpy.where(flength == 0, topo, -9999).argmax(), topo.shape
    )

    # While we still have points higher than 0 calculate the runoff cells from this high point
    while topo[highest_point] > 0:
        # Get cells connected to current cell
        cells, end_reason, next_cell = get_next_cell([highest_point], flength, trip)
        # We reverse the order of cells because the last one seen is the closest to the ocean
        x, y = numpy.array(cells[::-1]).T
        # We take the distance according to trip in each cell to the next cell
        # We calculate the cumlative sum along the flow path
        flow_lengths = numpy.cumsum(ocean_distances[x, y])

        if end_reason == "junction":
            # If we hit a junction then we need to adjust all the values because we are not starting at 0
            flength[x, y] = flength[next_cell] + flow_lengths
        elif end_reason == "ocean":
            # The outlet is an ocean
            flength[x, y] = flow_lengths

        else:
            raise AssertionError(
                f"Unexpected end reason {end_reason} for river starting at {highest_point}"
            )
        # Calculate the next highest point
        highest_point = numpy.unravel_index(
            numpy.where(flength == 0, topo, -9999).argmax(), topo.shape
        )

    # Add nan values for the ocean
    flength[omsk == 1] = numpy.nan
    return flength


def calculate_distbox(flength, trip):
    # Set the distance to ocean inside the ocean to 0
    tmp_flength = numpy.where(~numpy.isnan(flength), flength, 0)

    diffs = get_differences_with_neighbors(tmp_flength, ensure_gradient=True)

    xx, yy = numpy.ix_(numpy.arange(trip.shape[0]), numpy.arange(trip.shape[1]))
    # For trip values in the 1 -> range we calculate the difference between the cell and the nexy downstream cell
    assert numpy.nanmin(trip) > 0
    assert numpy.nanmax(trip) < 9
    distbox = (
        diffs[numpy.where(~numpy.isnan(trip), trip - 1, 0).astype(int), xx, yy] * 1000
    )
    return distbox


def calculate_trip_outflow_values(trip, outflow_points, basins, omsk, rlat):
    trip = trip.copy()
    trip[outflow_points[:, 0], outflow_points[:, 1]] = 9

    x, y = outflow_points.T

    for k in range(len(x)):
        i, j = x[k], y[k]
        adjacent_coords = get_adjacent_coords(i, j, trip)
        ip1, jp1 = numpy.max(adjacent_coords, axis=0)
        im1, jm1 = numpy.min(adjacent_coords, axis=0)

        adjacent_basins_values = get_adjacent_values(i, j, basins)
        basbas = [
            adjacent_basins_values[1, 1],
            adjacent_basins_values[2, 2],
            adjacent_basins_values[1, 2],
            adjacent_basins_values[2, 1],
        ]

        adjacent_trip_values = get_adjacent_values(i, j, trip)
        trptrp = [
            adjacent_trip_values[1, 1],
            adjacent_trip_values[2, 2],
            adjacent_trip_values[1, 2],
            adjacent_trip_values[2, 1],
        ]

        if ((basbas == basbas[0]).sum() == 4) & ((trptrp == trptrp[0]).sum() == 4):
            # Count how many ocean cells are in each of the corners
            ocean_in_box = [
                (
                    ~numpy.isnan(get_adjacent_values(im1, jm1, basins))
                ).sum(),  # Upper left
                # Upper right
                (~numpy.isnan(get_adjacent_values(im1, jp1, basins))).sum(),
                # Lower Right
                (~numpy.isnan(get_adjacent_values(ip1, jp1, basins))).sum(),
                (
                    ~numpy.isnan(get_adjacent_values(ip1, jm1, basins))
                ).sum(),  # Lower Left
            ]
            li = numpy.argmax(ocean_in_box)
            if li == 0:
                trip[i, j] = 9
                trip[ip1, j] = 7
                trip[ip1, jp1] = 8
                trip[i, jp1] = 1
            elif li == 1:
                trip[i, j] = 3
                trip[ip1, j] = 9
                trip[ip1, jp1] = 1
                trip[i, jp1] = 2
            elif li == 2:
                trip[i, j] = 4
                trip[ip1, j] = 5
                trip[ip1, jp1] = 9
                trip[i, jp1] = 3
            else:
                trip[i, j] = 5
                trip[ip1, j] = 6
                trip[ip1, jp1] = 7
                trip[i, jp1] = 9

    # We have modified the outflow points in the last part of code so outflow points are no longer correct
    # We need to find trip values equal to 9
    x, y = numpy.where(trip == 9)

    # We pad the arrays twice because we are looking at the 5 x 5 grid centered at the point where trip is equal to 9
    # If we don't pad twice then we will have problems on the borders
    padded_omsk = get_padded_array(omsk)  # pad 1 value
    padded_omsk = get_padded_array(padded_omsk)  # pad 2 value

    padded_trip = get_padded_array(trip)  # pad 1 value
    padded_trip = get_padded_array(padded_trip)  # pad 2 value

    for k in range(len(x)):
        i, j = x[k], y[k]
        # To get the five values centered on i, j we need to go from i - 2 -> i + 3
        # has we have padded twice we need to add 2 so the indexs match
        omsk_bx = padded_omsk[i - 2 + 2 : i + 3 + 2, j - 2 + 2 : j + 3 + 2]
        trip_bx = padded_trip[i - 2 + 2 : i + 3 + 2, j - 2 + 2 : j + 3 + 2]

        # If there is at least one ocean point then it is a coastal or river flow
        # Ocean values are 1 in omsk
        if numpy.sum(omsk_bx) > 0:
            # The first 200 basins are river flow else it is coastal flow
            if basins[i, j] < 200:
                trip[i, j] = 99
            else:
                trip[i, j] = 98
        elif numpy.sum(omsk_bx < 0.5) > 0:
            # This can only be an internal basin
            # Greenland is still a problem as we have coarse resolution coast lines
            # Thus anything north of 60deg N will be coastal flow
            if rlat[i, j] > 60:
                trip[i, j] = 98
            else:
                trip[i, j] = 97
        # Not sure we every hit this?
        else:
            raise AssertionError(
                "We have an ouflow point but we can not say if it is return flow or flow to the ocean"
            )
            print(
                "We have an ouflow point but we can not say if it is return flow or flow to the ocean"
            )
            print(i, j)
            print(omsk_bx)
            print(trip_bx)

    return trip


def calculate_dzz(topo, trip, distbox, omsk):

    height_differences = get_differences_with_neighbors(topo, ensure_gradient=True)

    xx, yy = numpy.ix_(numpy.arange(trip.shape[0]), numpy.arange(trip.shape[1]))
    # For trip values in the 1 -> range we calculate the difference between the cell and the nexy downstream cell
    dzz = height_differences[numpy.where(trip < 50, trip - 1, 0).astype(int), xx, yy]
    # For trip values corresponding to outlets we just take the topo value
    dzz = numpy.where(trip > 50, topo, dzz)
    # If values are less than 5 then clip them to "avoid unpleasant surprises"
    dzz = numpy.where(dzz > 5, dzz, 5)
    # Replace ocean values with 0
    # omsk:
    # - ocean = 1
    # - land = 0
    # Keep dzz values on land and everything else make 0
    dzz = numpy.where(omsk == 0, dzz, 0)
    return dzz


def calculate_topo_index(distbox, dzz, omsk):
    topoindex = numpy.sqrt(distbox ** 3.0 / (dzz * 10 ** 6))
    # omsk == 0 -> Land values
    # keep topoindex on land and NaN in the ocean
    topoindex = numpy.where(omsk == 0, topoindex, numpy.nan)
    return topoindex


def _set_attributes_ds_routing(ds_routing):
    ds_routing["nav_lon"].attrs = {
        "units": "degrees_east",
        "valid_min": -180.0,
        "valid_max": 180.0,
        "long_name": "Longitude",
    }
    ds_routing["nav_lat"].attrs = {
        "units": "degrees_north",
        "valid_min": -90.0,
        "valid_max": 90.0,
        "long_name": "Latitude",
    }
    ds_routing["trip"].attrs = {
        "axis": "TYX",
        "units": "",
        "long_name": "Direction of flow from each box",
        "associate": "nav_lat na",
    }
    ds_routing["basins"].attrs = {
        "axis": "TYX",
        "units": "",
        "long_name": "Basin number for each grid box",
        "associate": "nav_lat na",
    }
    ds_routing["topoind"].attrs = {
        "axis": "TYX",
        "units": "",
        "long_name": "Topographic index of the retention time",
        "associate": "nav_lat na",
    }
    ds_routing["hdiff"].attrs = {
        "axis": "TYX",
        "units": "m",
        "long_name": "Height difference with the outflow basin",
        "associate": "nav_lat na",
    }
    ds_routing["riverl"].attrs = {
        "axis": "TYX",
        "units": "m",
        "long_name": "River length to the next basin",
        "associate": "nav_lat na",
    }
    ds_routing["orog"].attrs = {
        "axis": "TYX",
        "units": "m",
        "long_name": "Orography",
        "associate": "nav_lat na",
    }
    ds_routing["disto"].attrs = {
        "axis": "TYX",
        "units": "km",
        "long_name": "Distance to ocean",
        "associate": "nav_lat na",
    }
    return ds_routing


def create_routing_netcdf(
    topo, trip, basins, topo_index, dzz, distbox, orog, flength, rlat, rlon
):
    ds_routing = xr.Dataset(
        coords={
            "nav_lon": (["y", "x"], rlon),
            "nav_lat": (["y", "x"], rlat),
        },
        data_vars={
            "trip": (["y", "x"], trip[::-1]),
            "basins": (["y", "x"], basins[::-1]),
            "topoind": (["y", "x"], topo_index[::-1]),
            "hdiff": (["y", "x"], dzz[::-1]),
            "riverl": (["y", "x"], distbox[::-1]),
            "orog": (["y", "x"], orog[::-1]),
            "disto": (["y", "x"], flength[::-1]),
            "topo": (["y", "x"], topo[::-1]),
        },
    )
    ds_routing = _set_attributes_ds_routing(ds_routing)
    return ds_routing


def create_bathy_paleo_orca(topo, custom_orca=None):
    # Transform topo to bathy and mask land
    bathy = topo * -1
    # Ocean values are now positive or equal to 0
    # Any strictly negative value is a continent
    # TODO we could pass omsk in for consitency
    bathy = numpy.where(bathy < 0, numpy.nan, bathy) # set land points to nans
    # Remap to nemo grid
    # Load grid
    ds_orca = custom_orca
    # Format default coordinate file or the one provided by user
    # If ds has no coords defined
    if not ds_orca.coords:
        listVar = list(ds_orca.data_vars)
        # finds lon and lat variable in ds_orca
        varLon = [var for var in listVar if "lon" in var.lower()]
        varLat = [var for var in listVar if "lat" in var.lower()]
        # Ensure there are only one lon and lat variables found and assign them as coords
        if len(varLon) == 1 and len(varLat) == 1:
            ds_orca = ds_orca.assign_coords(lon=(ds_orca[varLon[0]]))
            ds_orca = ds_orca.assign_coords(lat=(ds_orca[varLat[0]]))
        else:
            raise NameError(
                '''
                Can''t identify coordinates in "Coords" file. Several variables 
                containing "lon" and/or several variables containing "lat" found !
                '''
                )
    # Remove the missing value encoding fields as this causes problems because _FillValues is set to nan
    try:
        del ds_orca.lat.encoding["missing_value"]
    except KeyError:
        pass
    try:
        del ds_orca.lon.encoding["missing_value"]
    except KeyError:
        pass
    # create output dataset
    ds = xr.Dataset(
        coords={"lat": numpy.arange(-89.5, 90), "lon": numpy.arange(-179.5, 180)},
        data_vars={"Bathymetry": (["lat", "lon"], bathy)},
    )
    ds = ds.interp(ds_orca.coords, method="nearest", kwargs={"fill_value": None})
    ds.Bathymetry.values[ds.Bathymetry.values == 0] = 10 # set water points from value 0 to 10 in order to not be confused with land points that are nans and will get value 0 in next line
    ds.Bathymetry.values = numpy.where(
        numpy.isnan(ds.Bathymetry.values), 0, ds.Bathymetry.values
    ) # set land points from nans to value 0
    ds = ds.rename({"lat": "nav_lat", "lon": "nav_lon"})
    ds.nav_lat.attrs = {
        "standard_name": "latitude",
        "long_name": "latitude",
        "units": "degrees_north",
        "_CoordinateAxisType": "Lat",
    }
    ds.nav_lon.attrs = {
        "standard_name": "longitude",
        "long_name": "longitude",
        "units": "degrees_east",
        "_CoordinateAxisType": "Lon",
    }
    return ds


def create_topo_high_res(topo):
    # load new coordinates
    lon_vals = numpy.arange(-180 + (360.0 / 2160.0) / 2, 180, 360.0 / 2160)
    # Y values are funky so are stored in a file
    y_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "topo_high_res_y_vals.npy"
    )
    lat_vals = numpy.load(y_file)
    # create output dataset
    ds = xr.Dataset(
        coords={
            "latitude": numpy.arange(-89.5, 90),
            "longitude": numpy.arange(-179.5, 180),
        },
        data_vars={"RELIEF": (["latitude", "longitude"], topo)},
    )
    ds = ds.interp(
        {
            "latitude": lat_vals,
            "longitude": lon_vals,
        },
        method="linear",
    )

    ds.longitude.attrs = {
        "standard_name": "longitude",
        "long_name": "longitude",
        "units": "degrees",
        "axis": "X",
    }
    ds.latitude.attrs = {
        "standard_name": "latitude",
        "long_name": "latitude",
        "units": "degrees",
        "axis": "Y",
    }
    ds.RELIEF.attrs = {"long_name": "relief"}
    relief = ds.RELIEF.values.copy()
    # Ocean values are <= 0
    # < 0 chooses all ocean points not already set to 0
    relief[relief < 0] = 0
    relief[numpy.isnan(relief)] = 0
    ds.RELIEF.values = relief
    return ds


def create_soils(rlat, rlon, omsk):
    soil_color = numpy.ones(rlat.shape) * 4.0
    # Mask Oceans
    soil_color[omsk == 1] = 0
    soil_text = numpy.ones(rlat.shape) * 3.0
    # Mask Oceans
    soil_text[omsk == 1] = 0

    # create output dataset
    ds = xr.Dataset(
        coords={},
        data_vars={
            "nav_lon": (["y", "x"], rlon),
            "nav_lat": (["y", "x"], rlat),
            "soilcolor": (["y", "x"], soil_color[::-1]),
            "soiltext": (["y", "x"], soil_text[::-1]),
        },
    )
    ds["nav_lon"].attrs = {
        "units": "degrees_east",
        "valid_min": -180.0,
        "valid_max": 180.0,
        "long_name": "Longitude",
    }
    ds["nav_lat"].attrs = {
        "units": "degrees_north",
        "valid_min": -90.0,
        "valid_max": 90.0,
        "long_name": "Latitude",
    }
    ds["soilcolor"].attrs = {
        "axis": "YX",
        "units": "index",
        "long_name": "Soil color types from the Henderson-Sellers and Wilson dataset",
        "associate": "nav_lat nav_lon",
    }
    ds["soiltext"].attrs = {
        "axis": "YX",
        "units": "index",
        "long_name": "Soil texture from Zobler 86",
        "associate": "nav_lat nav_lon",
    }
    return ds


def run_routines(topo, latitudes, custom_orca=None):
    print(f" [c] {datetime.now()} Fixing Topo", flush=True)
    topo = fix_topo(topo, latitudes)
    print(f" [c] {datetime.now()} Calculating orogen", flush=True)
    orog = calculate_orog(topo)
    print(f" [c] {datetime.now()} Calculating Ocean Mask", flush=True)
    omsk = calculate_omsk(topo)
    print(f" [c] {datetime.now()} Calculating runoff directions", flush=True)
    trip = calculate_trip(topo, omsk)
    print(f" [c] {datetime.now()} Calculating curvilinear coords", flush=True)
    rlon, rlat = calculate_curvilinear_coordinates()
    print(f" [c] {datetime.now()} Calculating grid area", flush=True)
    area = calculate_area(rlat)
    print(f" [c] {datetime.now()} Calculating Runoff basins", flush=True)
    basins = calculate_basins(topo, trip, area)
    print(f" [c] {datetime.now()} Calculating distance to ocean", flush=True)
    ocean_distances = calculate_ocean_distances(trip, rlat)
    print(f" [c] {datetime.now()} Calculating length of rivers", flush=True)
    river_length = calculate_river_lengths(topo, trip, ocean_distances, omsk)
    print(f" [c] {datetime.now()} Calculating Distbox", flush=True)
    distbox = calculate_distbox(river_length, trip)
    print(f" [c] {datetime.now()} Calculating outflow points", flush=True)
    outflow_points = calculate_outflow_points(basins, trip)
    print(f" [c] {datetime.now()} Calculating Trip Values", flush=True)
    trip = calculate_trip_outflow_values(trip, outflow_points, basins, omsk, rlat)
    print(f" [c] {datetime.now()} Calculating dzz", flush=True)
    dzz = calculate_dzz(topo, trip, distbox, omsk)
    print(f" [c] {datetime.now()} Calculating topo_index", flush=True)
    topo_index = calculate_topo_index(distbox, dzz, omsk)

    print(f" [c] {datetime.now()} Creating routing file", flush=True)
    ds_routing = create_routing_netcdf(
        topo, trip, basins, topo_index, dzz, distbox, orog, river_length, rlat, rlon
    )
    print(f" [c] {datetime.now()} Creating Bathy file", flush=True)
    ds_bathy = create_bathy_paleo_orca(topo, custom_orca)
    print(f" [c] {datetime.now()} Creating Soils file", flush=True)
    ds_soils = create_soils(rlat, rlon, omsk)
    print(f" [c] {datetime.now()} Creating High Res file", flush=True)
    ds_topo_high_res = create_topo_high_res(topo)

    return ds_routing, ds_bathy, ds_soils, ds_topo_high_res
