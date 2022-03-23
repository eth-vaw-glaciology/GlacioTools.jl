# GlacioTools.jl

## Done already
General downloading routine `get_all_data()` (from `FourDAntarcticaSubglacialRouting.jl`) which can be used in three ways
```
get_all_data(myfile, destination_dir)
get_all_data(myfolder, destination_dir)
get_all_data(mydictionary, destination_dir)
```

## TODO
- download and read in Antarctica data (from `FourDAntarcticaSubglacialRouting.jl`)
- download and read data from Swiss alps (from `https://github.com/PTsolvers/FastIce.jl/tree/main/GeoData`)
- functions that help modifying raw data (e.g. smoothing), especially e.g. ice thickness

- documentation, in particular for downloads that require .netrc file
