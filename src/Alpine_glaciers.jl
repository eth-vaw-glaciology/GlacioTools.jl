"""
    fetch_glacier(name; datadir)

# Input
- `name` -- name of the glacier, used to name the files that are being saved
- `SGI_ID` -- unique ID of the glacier to find it in the dataset (found at https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/434697/00_TablesIllustrations%28updatedversion%29.pdf?sequence=39&isAllowed=y). [Note that in that dataset the separator is "/" (e.g. "N23/07") but now the separator "-" is used ("N23-07")]
- `datadir` -- path of the directory to store the downloaded data

# Output
- struct of type `GlacioTools.DataElevation` with fields ``x``, ``y``, ``z_bed``, ``z_surf`` and rotation matrix ``R``
"""
function fetch_glacier(name::String, SGI_ID::String; datadir::String)
    # Download data if needed
    fetch_data(datadir)
    # Extract galcier and crop
    geom_select(name, SGI_ID, datadir)
    # Save to non-GIS format
    extract_geodata(Float64, name, datadir)
    return load_elevation(joinpath(datadir,"alps/data_$(name).h5"))
end

"Move all files in one folder and remove files that are not needed."
function organise_folder(dir)
    files_to_move = [joinpath.(dir,"04_IceThickness_SwissAlps/IceThickness.tif")
                     joinpath.(dir,"08_SurfaceElevation_SwissAlps/",readdir(joinpath(dir,"08_SurfaceElevation_SwissAlps/")))
                     joinpath.(dir,"TLM_BB", filter(x->startswith(x,"swissTLM3D_TLM_BODEN") && (endswith(x,".shp") || endswith(x,".dbf")),readdir(joinpath(dir,"TLM_BB/"))))
                     joinpath.(dir,"TLM_BB", filter(x->startswith(x,"swissTLM3D_TLM_GLAMOS") && endswith(x,".dbf"),readdir(joinpath(dir,"TLM_BB/"))))]
    dirs_to_remove = joinpath.(dir,filter(x->isdir(joinpath(dir, x)),readdir(dir)))
    run(`mv $files_to_move $dir`)
    run(`rm -r $dirs_to_remove`)
    return
end

"""
    fetch_data(datadir::String)

Fetch the data from the internet

# Input
- `datadir` -- path of the directory to store the downloaded data
"""
function fetch_data(datadir::String)
    datas = Dict(   # Source: https://www.swisstopo.admin.ch/en/geodata/landscape/tlm3d.html#download
                    :swissTLM3D   => "https://data.geo.admin.ch/ch.swisstopo.swisstlm3d/swisstlm3d_2022-03/swisstlm3d_2022-03_2056_5728.shp.zip",
                    # Source: https://www.research-collection.ethz.ch/handle/20.500.11850/434697
                    :icethickness => "https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/434697/04_IceThickness_SwissAlps.zip?sequence=10&isAllowed=y",
                    :swissalti3D  => "https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/434697/08_SurfaceElevation_SwissAlps.zip?sequence=41&isAllowed=y"
                    )
    final_files = ["IceThickness.tif", "SwissALTI3D_r2019.tif", "swissTLM3D_TLM_GLAMOS.dbf",
                   "swissTLM3D_TLM_BODENBEDECKUNG_OST.shp", "swissTLM3D_TLM_BODENBEDECKUNG_WEST.shp",
                   "swissTLM3D_TLM_BODENBEDECKUNG_OST.dbf", "swissTLM3D_TLM_BODENBEDECKUNG_WEST.dbf"]
    # download
    sgi_dir = joinpath(datadir,"alps_sgi/")
    mkpath(sgi_dir)
    if any(.!isfile.(joinpath.(sgi_dir,final_files)))
        get_all_data(datas,sgi_dir)
        organise_folder(sgi_dir)
    end
    return
end

"""
    geom_select(name::String, SGI_ID::String, datadir; padding=100, do_save=true)

Select ice thickness, surface and bedrock elevation data for a given Alpine glacier based on SGI ID.

# Input
- `name::String`: input data file
- `SGI_ID::String`: desired data type for elevation data
- `datadir` -- path of the directory to store the downloaded data

# Optional keyword args
- `padding=100`: padding around the glacier geometry in meters
- `do_save=true`: save output to HDF5
"""
@views function geom_select(name::String, SGI_ID::String, datadir; padding=100, do_save=true)
    if isfile(joinpath(datadir,"alps/$(name)_IceThick_cr0.tif")) && isfile(joinpath(datadir,"alps/$(name)_SurfElev_cr.tif")) && isfile(joinpath(datadir,"alps/$(name)_BedElev_cr.tif"))
        return
    end

    # find glacier ID: there might be several polygons belonging to the same glacier
    df = DataFrame(DBFTables.Table(joinpath(datadir,"alps_sgi/swissTLM3D_TLM_GLAMOS.dbf")))
    ID = df[in([SGI_ID]).(df.SGI),:TLM_BODENB] # and not :UUID field!

    ##m3: alternative not needing DataFrames
    # dt = DBFTables.Table(joinpath(datadir,"alps_sgi/swissTLM3D_TLM_GLAMOS.dbf"))
    # inds = findall(dt.SGI .== SGI_ID)
    # IDs = dt.TLM_BODENB[inds]

    # read in global data
    print("Reading in global data... ")
    IceThick = read(Raster(joinpath(datadir,"alps_sgi/IceThickness.tif")))
    SurfElev = read(Raster(joinpath(datadir,"alps_sgi/SwissALTI3D_r2019.tif")))
    println("done.")

    count = 0
    IceThick_stack = []
    for id in ID
        count+=1
        # retrieve shape
        dftable = DataFrame(Shapefile.Table(joinpath(datadir,"alps_sgi/swissTLM3D_TLM_BODENBEDECKUNG_OST.shp")))
        if sum(in([id]).(dftable.UUID))==0
            dftable = DataFrame(Shapefile.Table(joinpath(datadir,"alps_sgi/swissTLM3D_TLM_BODENBEDECKUNG_WEST.shp")))
        end
        shape = dftable[in([id]).(dftable.UUID),:geometry][1] #m3: does shape need to be a Vector?
        # find ice thickness for polygon of interest (glacier), crop and add padding, using global data
        # push!(IceThick_stack, GlacioTools.crop_padded(IceThick, shape, padding)) #lr: does not correctly mask the "no ice" data.
        push!(IceThick_stack, mask_trim(IceThick, shape; pad=padding))
    end

    IceThick_cr = mosaic(first, IceThick_stack)

    # crop surface elevation to ice thickness data
    SurfElev_cr = Rasters.crop(SurfElev; to=IceThick_cr)

    # compute bedrock elevation
    IceThick_cr0 = replace_missing(IceThick_cr, 0.0)
    BedElev_cr   = SurfElev_cr .- IceThick_cr0

    print("Saving to file... ")
    # save
    if do_save
        if isdir(joinpath(datadir,"alps"))==false mkdir(joinpath(datadir,"alps")) end
        write(joinpath(datadir,"alps/$(name)_IceThick_cr0.tif"), IceThick_cr0)
        write(joinpath(datadir,"alps/$(name)_SurfElev_cr.tif") , SurfElev_cr)
        write(joinpath(datadir,"alps/$(name)_BedElev_cr.tif")  , BedElev_cr)
    end
    println("done.")

    return IceThick_cr0, SurfElev_cr, BedElev_cr
end

"""
    extract_geodata(type::DataType, dat_name::String)

Extract geadata and return bedrock and surface elevation maps, spatial coords and bounding-box rotation matrix.

# Inputs
- `type::DataType`: desired data type for elevation data
- `dat_name::String`: name of the glacier
"""
@views function extract_geodata(type::DataType, dat_name::String, datadir)
    if isfile(joinpath(datadir,"alps/data_$(dat_name).h5"))
        return
    end
    println("Starting geodata extraction ...")
    println("- load the data")
    file1     = joinpath(datadir,"alps/$(dat_name)_IceThick_cr0.tif")
    file2     = joinpath(datadir,"alps/$(dat_name)_BedElev_cr.tif"  )
    z_thick_i = reverse(Raster(file1),dims=2)
    z_bed_i   = reverse(Raster(file2),dims=2)
    xy        = DimPoints(dims(z_thick_i, (X, Y)))
    (x,y)     = (first.(xy), last.(xy))
    xmin,xmax = extrema(x)
    ymin,ymax = extrema(y)
    # center data in x,y plane
    x       .-= 0.5*(xmin + xmax)
    y       .-= 0.5*(ymin + ymax)
    # TODO: a step here could be rotation of the (x,y) plane using bounding box (rotating calipers)
    # define and apply masks
    z_thick   = z_thick_i.data[:,:]
    z_bed     = z_bed_i.data[:,:]
    my_mask   = 1 .- (z_bed.==z_bed_i.missingval)
    z_bed[my_mask.==0] .= mean(my_filter(z_bed,my_mask))
    z_bed    .= reshape(z_bed,size(my_mask))
    # ground data in z axis
    zmin =  minimum(z_bed)
    z_bed   .-= zmin
    # ice surface elevation and average between bed and ice
    z_surf    = z_bed .+ z_thick
    z_avg     = z_bed .+ convert(type,0.5).*z_thick
    println("- perform least square fit")
    αx, αy = lsq_fit(my_filter(x,my_mask),my_filter(y,my_mask),my_filter(z_avg,my_mask))
    # normal vector to the least-squares plane
    # rotation axis - cross product of normal vector and z-axis
    nv = [-αx  ,-αy   ,1.0]; nv ./= norm(nv)
    ax = [nv[2],-nv[1],0.0]; ax ./= norm(ax)
    # rotation matrix from rotation axis and angle
    R  = axis_angle_rotation_matrix(ax,acos(nv[3]))
    println("- save data to $(datadir)alps/data_$(dat_name).h5")
    h5open(joinpath(datadir,"alps/data_$(dat_name).h5"),"w") do io
        g = create_group(io, "glacier")
        g["x",compress=3]        = convert.(type,x.data)
        g["y",compress=3]        = convert.(type,y.data)
        g["z_bed",compress=3]    = convert.(type,z_bed)
        g["z_surf",compress=3]   = convert.(type,z_surf)
        g["x_offset",compress=3] = convert(type, 0.5*(xmin + xmax))
        g["y_offset",compress=3] = convert(type, 0.5*(ymin + ymax))
        g["z_offset",compress=3] = convert(type, zmin)
        g["R",compress=3]        = convert.(type,R)
    end
    println("done.")
    return
end
