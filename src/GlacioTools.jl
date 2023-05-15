module GlacioTools
using Downloads, Rasters, DBFTables, DataFrames, Shapefile,
    Statistics, Interpolations, LinearAlgebra, HDF5,
    GeoInterface, ArchGDAL, DelimitedFiles
using Printf

export get_all_data, fetch_Antarctica, fetch_glacier, fetch_data, geom_select, extract_geodata

# define some useful constants
const rhoi = 917
const rhow = 1000
const g    = 9.81
const day  = 24*60*60
const year = 365*day

include("downloading_helpers.jl")
include("Rasters_helpers.jl")
include("Antarctica.jl")
include("Alpine_glaciers.jl")
include("misc.jl")

end # module
