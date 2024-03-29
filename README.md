# GlacioTools.jl

A package to assist in downloading and reading data, particularly glaciology relevant data from Antarctica and a few Alpine glaciers.

[![Build Status](https://github.com/eth-vaw-glaciology/GlacioTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/eth-vaw-glaciology/GlacioTools.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/eth-vaw-glaciology/GlacioTools.jl/branch/main/graph/badge.svg?token=WDSZG7V5R3)](https://codecov.io/gh/eth-vaw-glaciology/GlacioTools.jl)


## General download
The general downloading routine `get_all_data()` is flexible in its use:
```
get_all_data(myfile, datadir)
get_all_data(myfolder, datadir)
get_all_data(myurl, datadir)
get_all_data(mydictionary, datadir)
```
The file is saved to the folder `datadir`. If the input is a file or folder path, the file is simply copied. The input can also be a dictionary, e.g.:
```
datas = Dict(:dataset1 => "https://ex1.mat",
             :dataset2 => ["https://ex2.zip", "https://ex3.zip"])
```

### Password protected download
If the download requires a password, as e.g. for the [Antarctica
Bedmachine
dataset](https://urs.earthdata.nasa.gov/oauth/authorize?client_id=_JLuwMHxb2xX6NwYTb4dRA&response_type=code&redirect_uri=https%3A%2F%2Fn5eil01u.ecs.nsidc.org%2FOPS%2Fredirect&state=aHR0cDovL241ZWlsMDF1LmVjcy5uc2lkYy5vcmcvTUVBU1VSRVMvTlNJREMtMDc1Ni4wMDIvMTk3MC4wMS4wMS9CZWRNYWNoaW5lQW50YXJjdGljYV8yMDIwLTA3LTE1X3YwMi5uYw),
the automatic download requires a `~/.netrc` file. In case of the bedmachine this file should contain a line
```
machine urs.earthdata.nasa.gov login myusername password 1234567
```
where `myusername` and `1234567` should be replaced by username and password, respectively. Note that no super-sensitive logins should be stored in such a way as the `.netrc` is not encrypted.

## Specific datasets

### Alpine glaciers
Download and read in the bed and surface topography of a selection of alpine glaciers (with routines from [FastIce.jl/GeoData](https://github.com/PTsolvers/FastIce.jl/tree/main/GeoData)).
```
SGI_ID = "B43-03"
datadir = joinpath(@__DIR__,"mydata/")
data = fetch_glacier("Rhone", SGI_ID; datadir)
```
where `data` is a struct with entries `x`,`y`,`z_bed`,`z_surf` and `R` (rotation matrix). The `SGI_ID` number identifies the glacier in the dataset and can be found [here](https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/434697/00_TablesIllustrations(updatedversion).pdf?sequence=39&isAllowed=y). This is an example dictionary containing the `SGI_ID`s for six selected glaciers:
```
SGI_IDS = Dict( "Rhone"       => "B43-03",
                "Aletsch"     => "B36-26",
                "PlaineMorte" => "A55f-03",
                "Morteratsch" => "E22-03",
                "Arolla"      => "B73-14",
                "ArollaHaut"  => "B73-12")
```

An alternative workflow is to first fetch the data for the entire Alps (the download step) and then select the geometry of a list of glacier given in the `SGI_IDS` dict automatically as following:
```
datadir = joinpath(@__DIR__,"mydata/")
fetch_data(datadir)
data = geom_select.(keys(SGI_IDS), values(SGI_IDS), datadir);
```

The data is downloaded from:
- [Swisstopo swissTLM3D](https://www.swisstopo.admin.ch/en/geodata/landscape/tlm3d.html#download)
- [Grab 2020, Swiss Glacier Thickness – Release 2020 (ETH Research Collection)](https://www.research-collection.ethz.ch/handle/20.500.11850/434697) (`04_IceThickness_SwissAlps.zip` and `08_SurfaceElevation_SwissAlps.zip`)

### Antarctica
Download and read Antarctica topography data
```
fetch_antarctica([:bedmachine]; datadir)
bm = read_bedmachine(datadir, thin=1)[1]
# etc.
```

## TODO
Helper functions:
- functions that help modifying raw data (e.g. smoothing), especially e.g. ice thickness

Readers:
...

Download:
- be clever on what to do with big datasets of which only a small part is used.
- be clever to not unzip everything if not everything is needed
  (swisstopo is such a case)

## Side notes

### Raster data as point or area data -- and -- as center coordinates

Raster data comes in two flavours, GDAL calls them Pixel-Is-Point or
Pixel-Is-Area.  To add to the confusion these two terms mix (1)
whether the data should be considered areal data or point data and (2)
how the coordinate system is defined (the former point to cell centres
whereas the latter to one of the edges of a pixel).  This can lead to
half a pixel offset.

In my opinion, madness lies in Pixel-Is-Area (not because of the area
interpretation but because of the half-cell coordinate offset) and any
data coming in as such should be converted to Pixel-Is-Point.

We need to make sure to be consistent there.

Rasters.jl and DimensionalData.jl can do either:
https://rafaqz.github.io/DimensionalData.jl/stable/api/#Sampling and
also disentangles the two concepts GDAL mixes.  We want all data
`DimensionalData.Dimensions.LookupArrays.Points` (which automatically
means that it is "center" coordinates).


Refs:
- https://www.ncei.noaa.gov/products/grid-cell-registration
- https://github.com/opengeospatial/coverage-implementation-schema/issues/6#issuecomment-768532086
- https://trac.osgeo.org/geotiff/wiki/RefiningGeoTIFF
- https://github.com/opengeospatial/ogcapi-coverages/issues/92
