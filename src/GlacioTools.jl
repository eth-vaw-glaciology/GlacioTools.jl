module GlacioTools
using Revise, Infiltrator # for development
using Downloads, Rasters, DBFTables, DataFrames, Shapefile
using Printf

export get_all_data, fetch_Antarctica, fetch_glacier

# define some useful constants -- or is this more annoying than useful?
const rhoi = 917
const rhow = 1000
const g    = 9.81
const day = 24*60*60
const year = 365*day

include("downloading_helpers.jl")
include("GeoData_helpers.jl")
include("Antarctica.jl")
include("Alpine_glaciers.jl")

end # module
