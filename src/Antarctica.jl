# FUNCTIONS TO DOWNLOAD AND READ IN ANTARCTICA DATA

using MAT

const box_antarctica = Box((-2750000, 2780000+1), (-2200000, 2300000+1))

# polar stereographic
const crs = EPSG(3031)
const crs_latlon = EPSG(4326)


vaw_url = "https://people.ee.ethz.ch/~werderm/4d-data-9xWArBUYVr/"
datas = Dict(:bedmachine => "https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0756.002/1970.01.01/BedMachineAntarctica_2020-07-15_v02.nc",   # requires password, i.e. .netrc file
             #:rema100 => "http://data.pgc.umn.edu/elev/dem/setsm/REMA/mosaic/v1.1/100m/REMA_100m_dem.tif",
             :basal_amery_2km => vaw_url * "goldberg/Amery_basal_melt/Amery_basal_melt_2km.mat",
             :basal_amery_5km => vaw_url * "goldberg/Amery_basal_melt/Amery_basal_melt_5km.mat",
             :bedmap2 => ["https://secure.antarctica.ac.uk/data/bedmap2/bedmap2_tiff.zip",
                          "https://secure.antarctica.ac.uk/data/bedmap2/bedmap2_readme.txt"],
             :lebrocq_flux => "ftp://ftp.quantarctica.npolar.no/Quantarctica3/Glaciology/Subglacial%20Water%20Flux/SubglacialWaterFlux_Modelled_1km.tif",
             :lakes_WrightSiegert => ["ftp://ftp.quantarctica.npolar.no/Quantarctica3/Glaciology/Subglacial%20Lakes/SubglacialLakes_WrightSiegert.shp",
                                      "ftp://ftp.quantarctica.npolar.no/Quantarctica3/Glaciology/Subglacial%20Lakes/SubglacialLakes_WrightSiegert.dbf"],
             :lakes_amery_hogg => vaw_url * "hogg/annas-lakes.tar.gz",
             :lakes_thwaites_malczyk => "https://4d-antarctica.org/wp-content/uploads/2021/01/Malczyk_etal_2020_data_v2.tar",
             :glads_amery_dow => [vaw_url * a for a in ["dow/gridded-outputs/Amery_ch_dis_hc_grid.csv",
                                                        "dow/gridded-outputs/Amery_ch_dis_lc_grid.csv",
                                                        "dow/gridded-outputs/Amery_X.csv",
                                                        "dow/gridded-outputs/Amery_Y.csv",
                                                        ]],
             :gls_measures => ["https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0709.002/1992.02.07/GroundingLine_Antarctica_v02.shp",
                               "https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0709.002/1992.02.07/GroundingLine_Antarctica_v02.dbf"],
             :basins_measures => ["https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0709.002/1992.02.07/Basins_Antarctica_v02.shp",
                                  "https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0709.002/1992.02.07/Basins_Antarctica_v02.dbf",
                                  "https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0709.002/1992.02.07/Basins_Antarctica_v02.shx",
                                  "https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0709.002/1992.02.07/Basins_Antarctica_v02.prj"],
             :basins_IMBIE => ["https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0709.002/1992.02.07/Basins_IMBIE_Antarctica_v02.shp",
                                  "https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0709.002/1992.02.07/Basins_IMBIE_Antarctica_v02.dbf",
                                  "https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0709.002/1992.02.07/Basins_IMBIE_Antarctica_v02.shx",
                                  "https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0709.002/1992.02.07/Basins_IMBIE_Antarctica_v02.prj"],
             :modis_moa2009 => ["https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0593_moa2009_v02/geotiff/moa125_2009_hp1_v02.0.tif.gz"]
             )


"""
    fetch_antarctica(keys=nothing; datadir="data/antarctica")

Download Antarctica data.  If only parts of the `GlacioTools.datas` should be downloaded,
then specify the keys in `datasets`.

The reading of the data can then be done with the appropriate file-reader.
"""
function fetch_antarctica(datasets=nothing; datadir="data/antarctica", kws...)
    dd = copy(datas)

    # keep only specified keys, default keep everything
    if !isnothing(datasets)
        for s in keys(dd) !in(s, datasets) ? delete!(dd, s) : nothing end
    end
    length(dd)==0 && return nothing
    mkpath(datadir)

    # download
    get_all_data(dd, datadir; kws...)

    return
end

"""
    read_bedmachine(datadir, thin=1; make_rmask=true, make_groundingline=true)

Reads the BedMachine Antractica dataset.

Also, by default, creates a routing mask and a mask of points just
land-wards of the routing mask.
"""
function read_bedmachine(datadir, thin=1; make_rmask=true, make_groundingline=true)
    nc = RasterStack(datadir * "/BedMachineAntarctica_2020-07-15_v02.nc")
    gas = []
    for k in [:bed, :errbed, :surface, :firn, :source, :mask]
        ga = nc[k]
        ga = (reverse(ga[1:thin:end, 1:thin:end], dims=2))[box_antarctica...] # this also loads it into memory
        if k==:errbed # this is {Missing, Int16}
            ga = replace_missing(ga, -9999)
            data = convert(Matrix{Float32}, ga.data)
            ga = replace_missing(Raster(data; ga.dims, ga.name, ga.refdims, metadata=Rasters.metadata(ga), missingval=-9999), NaN)
        end

        ga = if k==:bed || k==:surface || k==:firn
            # remove "missing" for bed and surface
            replace_missing(ga, NaN)
        elseif k==:source || k==:mask
            replace_missing(ga, -1)
        else
            ga
        end
        push!(gas, ga)

        if k==:mask && make_rmask
            # Add mask where (sub)glacial water should be routed/flowing:
            # grounded ice or lake Vostok or ice-free-land
            rmask = (ga.==2) .| (ga.==4) .| (ga.==1);
            push!(gas, Raster(rmask, name=:rmask))

            # make a mask which has all routing-points bordering on non-routing points
            # i.e. in a loose sense the grounding-line
            if make_groundingline
                groudingline = rmask .* false
                for IJ in CartesianIndices(rmask)
                    rmask[IJ]==false && continue # don't process water-points
                    for ij in iterate_D9(IJ, rmask)
                        ij==IJ && continue # don't process the point itself
                        if rmask[ij]==false # found water
                            groudingline[IJ] = true
                            continue
                        end
                    end
                end
                push!(gas, Raster(groudingline, name=:groundingline))
            end
        end
    end

    # make dims a range:
    x, y = dims(gas[1])

    return RasterStack(gas..., metadata=Rasters.metadata(nc), dims=(x, y)), nc
end

function read_bedmap2(datadir)
    fls = Dict(:bed => datadir*"/bedmap2_tiff/bedmap2_bed.tif",
               :surface => datadir*"/bedmap2_tiff/bedmap2_surface.tif")

    ga = Any[k=>reverse(Raster(v)[:,:,1],dims=2) for (k,v) in fls]

    errbed = reverse(Raster(datadir*"/bedmap2_tiff/bedmap2_grounded_bed_uncertainty.tif", missingval=typemax(UInt16))[:,:,1], dims=2)
    push!(ga, :errbed => errbed)
    rmask = Raster(datadir*"/bedmap2_tiff/bedmap2_icemask_grounded_and_shelves.tif")[:,end:-1:1,1].==0
    push!(ga, :rmask=>rmask)

    # make dims a range and also change their Sampling from Intervals to Points
    x, y = dims(errbed)
    dx = Int(step(x)); @assert step(y)==dx
    mod = x.val
    mod = Projected(;mod.order, mod.span, sampling=Rasters.DimensionalData.Dimensions.LookupArrays.Points(), mod.crs, mod.mappedcrs)
    D = (X(Int(x[1])+500:dx:round(Int,x[end])+500, mode=mod, metadata=metadata(x)),
         Y(Int(y[1])+500:dx:round(Int,y[end])+500, mode=mod, metadata=metadata(y)))
    @assert D[1].val .- 500==round.(Int, x.val) && D[2].val .- 500 == round.(Int, y.val)

    return RasterStack((;ga...), dims=D)[box_antarctica...]
end

## Quantartica datasets
#######################

read_lebrocq_flux(datadir) = Raster(datadir * "/SubglacialWaterFlux_Modelled_1km.tif")[1:end,1:end,1]
# function read_lebrocq_flux()
#     out = GDALarray(datadir * "/SubglacialWaterFlux_Modelled_1km.tif")[1:end,1:end,1]
#     out[out.==-9999] .= NaN
#     return Raster(out.data, dims=(X(dims(out)[1].val), Y(dims(out)[2].val)), name=:lebrocq)
# end

function read_wrightsiegert_lakes(datadir)
    # these are points
    geoms = Shapefile.shapes(Shapefile.Table(datadir * "/SubglacialLakes_WrightSiegert.shp"))
    return hcat([round.(Int,i) for i in GeoInterface.coordinates.(geoms)]...)
end

function read_gl_measures(datadir)
    # these are points
    geoms = Shapefile.shapes(Shapefile.Table(datadir * "/GroundingLine_Antarctica_v02.shp"))
    return hcat([round.(Int,i) for i in GeoInterface.coordinates.(geoms)[1][639][1]]...)
end

function read_basins_measures(datadir)
    # these are points
    geoms = Shapefile.shapes(Shapefile.Table(datadir * "/Basins_Antarctica_v02.shp"))
    return hcat([round.(Int,i) for i in GeoInterface.coordinates.(geoms)[1][639][1]]...)
end



## 4D Antarctica datasets
#########################

"""
    read_basal_melt_amery(D, datadir)

The `D` are the `dims` of an Antarctic dataset
"""
function read_basal_melt_amery(D, datadir)
    fl = joinpath(datadir, "Amery_basal_melt_2km.mat")
    handle = matopen(fl)
    dat = read(handle, "Amery_basal_melt_data")
    close(handle)

    mr = convert(Matrix{Float32}, dat["Melt_rate"]')
    mr[isnan.(mr)] .= 0
    mr[mr.<=0] .= 0 # TODO: are there freeze-on regions?

    # make into a Raster
    xx, yy = D
    x,y = dat["X_coord"][1,:], dat["Y_coord"][:,1]
    dx = Int(x[2]-x[1])
    xx = X(Int(x[1]):dx:Int(x[end]), mode=xx.val, metadata=metadata(xx))
    yy = Y(Int(y[1]):dx:Int(y[end]), mode=yy.val, metadata=metadata(yy))
    return Raster(mr/year, name=:melt, dims=(xx,yy)) # convert to m/s
end

function read_hogg_lakes(datadir)
    lakes_xy = Dict()
    for fl in readdir(joinpath(datadir,"annas-lakes/shapefiles"))
        if startswith(fl, "amery_lake") && splitext(fl)[2]==".shp"
            geoms = Shapefile.shapes(Shapefile.Table(datadir * "/annas-lakes/shapefiles/$fl"))[1]
            xy = hcat(GeoInterface.coordinates(geoms)[1][1]...)
            x,y = xy[1,:], xy[2,:]
            z = zeros(Float64, length(geoms.points))
            # transform to polar stereo (note that geoms cannot be directly transformed)
            ArchGDAL.createcoordtrans(crs_latlon, crs) do transform
                ArchGDAL.transform!(x, y, z, transform)
            end
            lakes_xy[Symbol(split(splitext(fl)[1], '_')[end])] = (x,y)
        end
    end
    return lakes_xy
end

function read_malczyk_lakes(datadir)
    lakes_masks = []
    for fl in readdir(datadir)
        if startswith(fl, "Thw") && splitext(fl)[2]==".tif"
            push!(lakes_masks, Raster(joinpath(datadir, fl))[1:end,1:end,1])
        end
    end
    out = deepcopy(lakes_masks[1])
    for lm in lakes_masks[2:end]
        out[lm.==1] .= 1
    end
    out[out.==0] .= NaN
    return out, lakes_masks
end

function read_glads_dow(datadir)
    out = []
    fls = [splitdir(f)[2] for f in datas[:glads_amery_dow]]
    # get XY
    x = convert(Matrix{Float32}, readdlm(joinpath(datadir, "Amery_X.csv"), ','))[1,:]
    y = convert(Matrix{Float32}, readdlm(joinpath(datadir, "Amery_Y.csv"), ','))[:,1]
    dx = Int(x[2]-x[1])
    xx = X(Int(x[1]):dx:Int(x[end]))
    yy = Y(Int(y[1]):dx:Int(y[end]))
    # get values
    dat = convert(Matrix{Float32}, readdlm(joinpath(datadir, "Amery_ch_dis_hc_grid.csv"), ',')')
    glads_hc = Raster(dat, name=:glads_hc, dims=(xx,yy))
    dat = convert(Matrix{Float32}, readdlm(joinpath(datadir, "Amery_ch_dis_lc_grid.csv"), ',')')
    glads_lc = Raster(dat, name=:glads_lc, dims=(xx,yy)) # convert to m/s

    return (glads_lc=glads_lc, glads_hc=glads_hc)
end
