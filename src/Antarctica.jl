# FUNCTIONS TO DOWNLOAD AND READ IN ANTARCTICA DATA

const box_antarctica = Box((-2750000, 2780000+1), (-2200000, 2300000+1))

# Note: To download the data the ISG shares need to be mounted (as described here https://vawiki.ethz.ch/vaw/informatics:samba_for_linux?s[]=samba
function fetch_Antarctica(spec=nothing;
                                destination_dir="data/Antarctica/",
                                bedmachine_thin=1)
    vaw_url = "https://people.ee.ethz.ch/~werderm/4d-data-9xWArBUYVr/"
    datas = Dict(:bedmachine => "https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0756.002/1970.01.01/BedMachineAntarctica_2020-07-15_v02.nc",   # requires password, i.e. .netrc file
                 #:bedmachine => "file://" * homedir() * "/itet-stor/glazioarch/GlacioData/BedMachine_Antarctica/168596330/BedMachineAntarctica_2019-11-05_v01.nc", # requires access to IGS share
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
                                   "https://n5eil01u.ecs.nsidc.org/MEASURES/NSIDC-0709.002/1992.02.07/GroundingLine_Antarctica_v02.dbf"]
                 )

    # keep only specified keys, default keep everything
    if !isnothing(spec)
        for s in keys(datas) !in(s, spec) ? delete!(datas, s) : nothing end
    end

    # download
    get_all_data(datas, destination_dir)

    # reading
    output = Dict()
    if haskey(datas, :bedmachine)
        topo, nc = read_bedmachine(destination_dir, bedmachine_thin)
        output[:bedmachine] = Dict(:topo => topo, :nc => nc)
    end
    if haskey(datas, :bedmap2)
        topo = read_bedmap2(destination_dir)
        output[:bedmap2] = topo
    end
    return output
end

function read_bedmachine(datadir, thin=1)
    nc = NCDstack(datadir * "/BedMachineAntarctica_2020-07-15_v02.nc") # this doesn't do anything: , childkwargs=(crs=crs,))
    # read arrays into memory and thin them, if desired
    gas = []
    for k in keys(nc)
        if k in [:mapping, :thickness,  :source, :geoid, :mask]
            continue # drop these fields
        end
        ga = (reverse(nc[k][1:thin:end, 1:thin:end], dims=2))[box_antarctica...] # this also loads it into memory
        ga = if k==:bed || k==:errbed || k==:surface || k==:firn
            # remove "missing" for bed, errbed and surface
            missing2nan(ga)
        else
            ga
        end
        push!(gas, ga)
    end
    mask = reverse(nc[:mask][1:thin:end, 1:thin:end], dims=2)

    # add mask where water is routed
    # route water over(under): grounded ice or lake Vostok or ice-free-land
    rmask = (mask.==2) .| (mask.==4) .| (mask.==1);
    push!(gas, Raster(rmask, name=:rmask)[box_antarctica...])

    # make dims a range:
    x, y = dims(gas[1])

    #dx = x[2]-x[1]; @assert y[2]-y[1]==dx
    #D = (X(x[1]:dx:x[end], mode=mode(x), metadata=metadata(x)),
    #     Y(y[1]:dx:y[end], mode=mode(y), metadata=metadata(y)))
    #@assert val(D[1])==val(x) && val(D[2])==val(y)
    # the double-RasterStack is necessary as the first cannot use the dims-kw:
    return RasterStack(RasterStack(gas..., metadata=nc.metadata), dims=(x, y)), nc
#    return RasterStack(gas..., metadata=nc.metadata), nc
end

function read_bedmap2(datadir)
    fls = Dict(:bed => datadir*"/bedmap2_tiff/bedmap2_bed.tif",
               :surface => datadir*"/bedmap2_tiff/bedmap2_surface.tif")

    ga = Any[k=>missing2nan(Raster.(v)[:,end:-1:1,1]) for (k,v) in fls]

    errbed = Raster(datadir*"/bedmap2_tiff/bedmap2_grounded_bed_uncertainty.tif", missingval=typemax(UInt16))[:,end:-1:1,1]
    push!(ga, :errbed => missing2nan(errbed))
    rmask = Raster(datadir*"/bedmap2_tiff/bedmap2_icemask_grounded_and_shelves.tif")[:,end:-1:1,1].==0
    push!(ga, :rmask=>rmask)

    # make dims a range and also change their Sampling from Intervals to Points
    x, y = dims(errbed)

    dx = Int(step(x)); @assert step(y)==dx
    mod = mode(x)
    mod = Projected(;mod.order, mod.span, sampling=Points(), mod.crs, mod.mappedcrs)
    D = (X(Int(x[1])+500:dx:round(Int,x[end])+500, mode=mod, metadata=metadata(x)),
         Y(Int(y[1])+500:dx:round(Int,y[end])+500, mode=mod, metadata=metadata(y)))
    @assert val(D[1]).-500==round.(Int,val(x)) && val(D[2]).-500==round.(Int,val(y))

    RasterStack(RasterStack((;ga...)), dims=D)[box_antarctica...]
end

# function read_basal_melt_amery(D)
#     fl = joinpath(datadir, "Amery_basal_melt_2km.mat")
#     handle = matopen(fl)
#     dat = read(handle, "Amery_basal_melt_data")
#     close(handle)

#     mr = convert(Matrix{Float32}, dat["Melt_rate"]')
#     mr[isnan.(mr)] .= 0
#     mr[mr.<=0] .= 0 # TODO: are there freeze-on regions?

#     # make into a Raster
#     xx, yy = D
#     x,y = dat["X_coord"][1,:], dat["Y_coord"][:,1]
#     dx = Int(x[2]-x[1])
#     xx = X(Int(x[1]):dx:Int(x[end]), mode=mode(xx), metadata=metadata(xx))
#     yy = Y(Int(y[1]):dx:Int(y[end]), mode=mode(yy), metadata=metadata(yy))
#     return Raster(mr/year, name=:melt, dims=(xx,yy)) # convert to m/s
# end

# function read_lebrocq_flux()
#     out = GDALarray(datadir * "/SubglacialWaterFlux_Modelled_1km.tif")[1:end,1:end,1]
#     out[out.==-9999] .= NaN
#     return Raster(out.data, dims=(X(dims(out)[1].val), Y(dims(out)[2].val)), name=:lebrocq)
# end

# function read_WrightSiegert_lakes()
#     # these are points
#     geoms = Shapefile.shapes(Shapefile.Table(datadir * "/SubglacialLakes_WrightSiegert.shp"))
#     return [round.(Int,i) for i in coordinates.(geoms)]
# end

# function read_gl_measures()
#     # these are points
#     geoms = Shapefile.shapes(Shapefile.Table(datadir * "/GroundingLine_Antarctica_v02.shp"))
#     return hcat([round.(Int,i) for i in coordinates.(geoms)[1][639][1]]...)
# end

# function read_hogg_lakes()
#     lakes_xy = Dict()
#     for fl in readdir(joinpath(datadir,"annas-lakes/shapefiles"))
#         if startswith(fl, "amery_lake") && splitext(fl)[2]==".shp"
#             geoms = Shapefile.shapes(Shapefile.Table(datadir * "/annas-lakes/shapefiles/$fl"))[1]
#             xy = hcat(coordinates(geoms)[1][1]...)
#             x,y = xy[1,:], xy[2,:]
#             z = zeros(Float64, length(geoms.points))
#             # transform to polar stereo (note that geoms cannot be directly transformed)
#             ArchGDAL.createcoordtrans(crs_latlon, crs2) do transform
#                 ArchGDAL.transform!(x, y, z, transform)
#             end
#             lakes_xy[Symbol(split(splitext(fl)[1], '_')[end])] = (x,y)
#         end
#     end
#     return lakes_xy
# end

# function read_malczyk_lakes()
#     lakes_masks = []
#     for fl in readdir(datadir)
#         if startswith(fl, "Thw") && splitext(fl)[2]==".tif"
#             push!(lakes_masks, GDALarray(joinpath(datadir, fl))[1:end,1:end,1])
#         end
#     end
#     out = deepcopy(lakes_masks[1])
#     for lm in lakes_masks[2:end]
#         out[lm.==1] .= 1
#     end
#     out[out.==0] .= NaN
#     return out, lakes_masks
# end

# function read_glads_dow()
#     out = []
#     fls = [splitdir(f)[2] for f in datas[:glads_amery_dow]]
#     # get XY
#     x = convert(Matrix{Float32}, readdlm(joinpath(datadir, "Amery_X.csv"), ','))[1,:]
#     y = convert(Matrix{Float32}, readdlm(joinpath(datadir, "Amery_Y.csv"), ','))[:,1]
#     dx = Int(x[2]-x[1])
#     xx = X(Int(x[1]):dx:Int(x[end]))
#     yy = Y(Int(y[1]):dx:Int(y[end]))
#     # get values
#     dat = convert(Matrix{Float32}, readdlm(joinpath(datadir, "Amery_ch_dis_hc_grid.csv"), ',')')
#     glads_hc = Raster(dat, name=:glads_hc, dims=(xx,yy))
#     dat = convert(Matrix{Float32}, readdlm(joinpath(datadir, "Amery_ch_dis_lc_grid.csv"), ',')')
#     glads_lc = Raster(dat, name=:glads_lc, dims=(xx,yy)) # convert to m/s

#     return (glads_lc=glads_lc, glads_hc=glads_hc)
# end
