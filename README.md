# GlacioTools.jl

A package to assist in downloading and reading data, particularly glaciology relevant data from Antarctica and a few Alpine glaciers.

## General download
The general downloading routine `get_all_data()` (from `FourDAntarcticaSubglacialRouting.jl`) is flexible in its use:
```
get_all_data(myfile, destination_dir)
get_all_data(myfolder, destination_dir)
get_all_data(myurl, destination_dir)
get_all_data(mydictionary, destination_dir)
```
The file is saved to the folder `destination_dir`. If the input is a file or folder path, the file is simply copied. The input can also be a dictionary, e.g.:
```
datas = Dict(:example1 => "https://ex1.mat",
             :example2 => "https://ex2.zip")
```

## Alpine glaciers
Download and read in the bed and surface topography of a selection of alpine glaciers (with routines from `FastIce.jl/GeoData`).
```
data = fetch_glacier("Rhone"; destination_dir)
```
where `data` is a struct with entries `x`,`y`,`z_bed`,`z_surf` and `R` (rotation matrix). Possible glaciers are:

- "Rhone"
- "Aletsch"
- "PlaineMorte"
- "Morteratsch"
- "Arolla"
- "ArollaHaut"



## Antarctica
Download and read Antarctica topography data, as in `FourDAntarcticaSubglacialRouting.jl`.
```
data = fetch_Antarctica([:bedmachine]; destination_dir)
```
So far only the :bedmachine is implemented.

## TODO
- make the reading functions for other Antarctica data work
- functions that help modifying raw data (e.g. smoothing), especially e.g. ice thickness
- documentation for downloads that require .netrc file
