import numpy
import xarray as xr


def load_pft(ds):
    # This function checks if data for PFTs has already been saved in db. 
    # IF NOT: generates a default configuration
    # IF YES: loads data in variables provided into html to render PFT table 

    if type(ds) == type(None): # If dataset ds doesn't exist
        percentVegetArrayClean = numpy.array([[0, 75, 25,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                                                [0, 15, 55,  0,  0,  0,  0,0,  0, 30,  0,  0,  0],
                                                [0,  0,  0,  0, 70, 30,  0,0,  0,  0,  0,  0,  0],
                                                [0,  0,  0,  0, 40, 30,  0,30,  0,  0,  0,  0,  0],
                                                [0,  0,  0,  0,  0, 30,  0,  0, 40, 30,  0,  0,  0]])
        boundsClean = numpy.array([15, 35, 50, 80, 90])

    else: # If dataset ds exist
        lenPFT = len(ds["veget"])
        lenLat = int(len(ds["lat"])/2)

        # From nc file, the 2 important data to retrieve are percent values for the different
        # PFTs and their locations that we define by extracting latitudinal bounds.
        
        percentVegetArray = numpy.zeros(shape=(18, lenPFT))  # 18 being max of bounds we can create in PFT step
        boundsVegetArray = numpy.zeros(shape=(18, lenPFT))
        percentVegetArray[:] = numpy.nan
        boundsVegetArray[:] = numpy.nan

        # Retrieve informations on percentages and bounds from netcdf pft file 
        # and store them in percentVegetArray and boundsVegetArray
        for indPFT in numpy.arange(0, lenPFT):
            percentVegetVect = numpy.zeros(shape=(lenLat))
            for j in numpy.arange(0, lenLat):
                percentVegetRaw = numpy.unique(ds["maxvegetfrac"][0, indPFT, j, :])
                percentVegetVect[j] = percentVegetRaw[~numpy.isnan(
                    percentVegetRaw)]  # get rid of nans
                if j == 1:
                    percentVeget = numpy.array([percentVegetVect[j]])
                    boundsVeget = numpy.array([90])
                elif j > 0 and percentVegetVect[j] != percentVegetVect[j-1]:
                    percentVeget = numpy.append(percentVeget, percentVegetVect[j])
                    boundsVeget = numpy.append(boundsVeget, ds["lat"][j]+0.25)
            for indBound in numpy.arange(0, len(boundsVeget)):
                percentVegetArray[indBound, indPFT] = percentVeget[indBound]
                boundsVegetArray[indBound, indPFT] = boundsVeget[indBound]

        # Create variable boundsClean containing the bounds and that will be provided in pft.html
        boundsClean = numpy.unique(boundsVegetArray)
        boundsClean = boundsClean[~numpy.isnan(boundsClean)].astype(int)

        lenBounds = len(boundsClean)

        # Truncate arrays to remove unecessary parts
        percentVegetArray = percentVegetArray[0:lenBounds, :]
        boundsVegetArray = boundsVegetArray[0:lenBounds, :]

        # Initialize final array percentVegetArrayClean
        percentVegetArrayClean = numpy.zeros(shape=(lenBounds, lenPFT))
        percentVegetArrayClean[:] = numpy.nan

        # From boundsVegetArray, percentVegetArray and boundsClean, reconstruct table of percentages 
        # percentVegetArrayClean that will be provided in pft.html
        for j in numpy.arange(0, lenPFT):
            for i in numpy.arange(0, lenBounds):
                if ~numpy.isnan(boundsVegetArray[i, j]):
                    ind = int(numpy.where(numpy.flip(boundsClean) == boundsVegetArray[i, j])[0])
                    percentVegetArrayClean[ind, j] = percentVegetArray[i, j]
            for ii in numpy.arange(0, lenBounds):
                if numpy.isnan(percentVegetArrayClean[ii, j]):
                    percentVegetArrayClean[ii, j] = percentVegetArrayClean[ii-1, j]

        percentVegetArrayClean = numpy.flipud(percentVegetArrayClean*100).astype(int)
    
    return boundsClean, percentVegetArrayClean


def generate_pft_netcdf(topo, latitudes, pft_values):
    # Attribute new PFTs. See lookup table to convert the 10 megabiomes of BIOME4 (Herold
    # GMD 2014, Harrison and Prentice Global Change Biology 2003) into the corresponding
    # 13 PFTs of ORCHIDEE (excluding n°13 and 14)
    #
    # The conversion is based on a rough comparison between the locations of the 10
    # megabiomes of a preindustrial BIOME4 simulation (see Herold GMD 2014) and the locations
    # of the 13 PFTs of the PFTmap_IPCC_1850.nc dataset used in preindustrial simulations
    # in the IPSL model. The same 10 megabiomes are then used by Herold GMD 2014 to create
    # a 55Ma global vegetation reconstruction.
    #
    # See also Lunt et al. GMDD 2016
    #
    # Abbreviation of ORCHIDEE PFTs:
    # BG     Bare ground  (k=1)
    # TBLE   Tropical broadleaved evergreen
    # TBLR   Tropical broadleaved raingreen
    # TNLE   Temperate needleleaf evergreen
    # TBLE2  Temperate broadleaved evergreen
    # TBLS   Temperate broadleaved summergreen
    # BNLE   Boreal needleleaf evergreen
    # BBLS   Boreal broadleaved summergreen
    # BNLS   Boreal needleleaf summergreen
    # C3     C3 grass
    # C4     C4 grass
    #
    # see also IPSL/modeles/ORCHIDEE/src_parameters/constantes_mtc.f90
    assert topo.shape == (360, 720)
    nb_lat = 360
    nb_lon = 720
    arr = numpy.zeros((13, nb_lat, nb_lon))

    # Make sure pft values are percentages
    assert numpy.max(pft_values) <= 100
    assert numpy.max(pft_values) > 1
    assert numpy.min(pft_values) >= 0

    # Convert percentages to 0 -> 1 range
    pft_values = pft_values / 100.0

    def _get_lat_index(val):
        return int((90 - val) / (180 / nb_lat))

    _from = 0
    for i in range(len(latitudes)):
        _to = latitudes[i]
        for pft_layer in numpy.argwhere(pft_values[i] > 0):
            # Set one hemisphere
            arr[pft_layer, _get_lat_index(_to) : _get_lat_index(_from), :] = pft_values[
                i, pft_layer
            ]
            # Set other hemisphere
            arr[
                pft_layer,
                nb_lat - _get_lat_index(_from) : nb_lat - _get_lat_index(_to),
                :,
            ] = pft_values[i, pft_layer]
        # Go to next latitude
        _from = _to

    # Convert oceans values to NaN values
    xx, yy = numpy.where(topo <= 0)
    arr[:, xx, yy] = numpy.nan

    ds = xr.Dataset(
        coords={
            "lat": numpy.arange(89.75, -90, -0.5),
            "lon": numpy.arange(-179.75, 180, 0.5),
            "time_counter": [1],
            "veget": numpy.arange(1, 14, dtype=numpy.int32),
        },
        data_vars={
            "maxvegetfrac": (
                ["time_counter", "veget", "lat", "lon"],
                arr[numpy.newaxis, ...],
            )
        },
    )
    ds.lat.attrs = {
        "bounds": "bounds_lat",
        "valid_max": 90.0,
        "long_name": "Latitude",
        "valid_min": -90.0,
        "units": "degrees_north",
        "axis": "Y",
    }
    ds.lon.attrs = {
        "bounds": "bounds_lon",
        "valid_max": 180.0,
        "long_name": "Longitude",
        "valid_min": -180.0,
        "axis": "X",
        "units": "degrees_east",
        "modulo": 360.0,
        "topology": "circular",
    }
    ds.time_counter.attrs = {
        "units": "years since 0-1-1",
        "calendar": "gregorian",
        "axis": "T",
    }
    ds.veget.attrs = {
        "valid_max": 13,
        "long_name": "Vegetation Classes",
        "valid_min": 1,
        "units": "-",
    }

    return ds
