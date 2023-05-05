# From FourDAntarcticaSubglacialRouting.jl #
############################################
get_x(ra::Union{Raster,RasterStack})  = dims(ra, X).val
get_y(ra::Union{Raster,RasterStack})  = dims(ra, Y).val
get_dx(ra::Union{Raster,RasterStack}) = diff(get_x(ra)[1:2])[1]
get_dy(ra::Union{Raster,RasterStack}) = diff(get_y(ra)[1:2][1:2])[1]

"""
    Box(x::Tuple, y::Tuple)

Get an index-tuple that can be used to crop a Raster to specified x and y coordinate ranges.

    ra[box...]

It is really just this:

    (X(x[1]..x[2]), Y(y[1]..y[2]))

TODO:
- probably not the right code.  Check GeoInterface.Extent
"""
Box(x::Tuple, y::Tuple) = (X(x[1]..x[2]), Y(y[1]..y[2]))

## Geography
function inbox(pt, box)
    xr,yr = Rasters.val.(box)
    (xr.left<=pt[1]<=xr.right) && (yr.left<=pt[2]<=yr.right)
end
function center(box)
    xr,yr = Rasters.val.(box)
    (xr.left+xr.right)/2, (yr.left+yr.right)/2
end

"""
    box2cartesian(box, ra::Raster)

Translates a coordinate box into CartesianIndices.
"""
box2cartesian(box, ra::Raster) = CartesianIndices(Rasters.DimensionalData.Dimensions.dims2indices(ra, box))

"""
    fix_GDAL_reversing_dims(raster; crs=dims(raster)[1].val.crs )

When an raster is processed by GDAL, it seems to like to switch the order of the Y-dim
and potentially add a Band-dim.  This tries to fix this.

Optionally pass in a `crs`
"""
function fix_GDAL_reversing_dims(raster; crs=dims(raster)[1].val.crs )
    # if a dim is reversed then re-reverse it
    # https://github.com/rafaqz/Rasters.jl/issues/385
    x,y = dims(raster)
    if y.val.order==DimensionalData.Dimensions.LookupArrays.ReverseOrdered()
        y = Y(LinRange(y[end],y[1], length(y)))
    else
        y = Y(LinRange(y[1],y[end], length(y)))
    end
    x = X(LinRange(x[1],x[end], length(x)))
    return Raster(reverse(raster.data[:,:],dims=2), (x,y); crs=crs, raster.name, raster.refdims, metadata=Rasters.metadata(raster))
end

function our_resample(raster; method=:bilinear, to, crs=dims(raster)[1].val.crs )
    out = resample(raster; to, method)
    return fix_GDAL_reversing_dims(out; crs)
    # if a dim is reversed then re-reverse it
    # https://github.com/rafaqz/Rasters.jl/issues/385
    # x,y = dims(out)
    # if y.val.order==DimensionalData.Dimensions.LookupArrays.ReverseOrdered()
    #     y = Y(LinRange(y[end],y[1], length(y)))
    # else
    #     y = Y(LinRange(y[1],y[end], length(y)))
    # end
    # x = X(LinRange(x[1],x[end], length(x)))
    # out = Raster(reverse(out.data,dims=2), (x,y); crs=crs, raster.name, raster.refdims, metadata=Rasters.metadata(raster))
    # return out
end

# From FastIce.jl/GeoData #
###########################

"""
    crop_padded(ra, po, pad=0)

Helper function to crop a raster `ra` to fully contain a polygon `po` plus some padding `pad`
(units of length).  This only works for 2D (3D is ignored).
"""
function crop_padded(ra, po, pad=0)
    ext = GeoInterface.extent(po)
    # extend Extent by pad, drop Z-dim
    ext = GeoInterface.Extent( (X=(ext.X[1]-pad, ext.X[2]+pad),
                                Y=(ext.Y[1]-pad, ext.Y[2]+pad)) )
    # crop
    return crop(ra, to=ext)
end

"""
    mask_trim(ra, po; pad=10)

Helper function to mask and trim a raster `ra` to fully contain a polygon `po` plus some
padding `pad` (units of length).  This only works for 2D (3D is ignored).
"""
mask_trim(ra, po; pad=10) = trim(mask(ra; with=po); pad=pad)

"Filter out all values of `A` based on `mask`."
my_filter(A, mask) = A[mask .!= 0]

"""
    lsq_fit(mask, zavg, xv2, yv2)

Linear least-squares fit of mean bedrock and surface data.
"""
function lsq_fit(x,y,z)
    # prepare input for least-squares regression (lsq)
    A       =  ones(length(x),3)
    B       = zeros(length(x),1)
    A[:,1] .= x; A[:,2] .= y; B .= z
    # lsq solve
    return (A'*A)\(A'*B)
end

"""
    axis_angle_rotation_matrix(ax, θ)

Get rotation matrix from rotation axis `ax` and angle `θ`.
"""
function axis_angle_rotation_matrix(ax, θ)
    return [cos(θ)+ax[1]^2*(1-cos(θ)) ax[1]*ax[2]*(1-cos(θ))       ax[2]*sin(θ)
            ax[2]*ax[1]*(1-cos(θ))    cos(θ) + ax[2]^2*(1-cos(θ)) -ax[1]*sin(θ)
           -ax[2]*sin(θ)              ax[1]*sin(θ)                       cos(θ)]
end

# From FastIce.jl/scripts3D/helpers3D_v5.jl #
#############################################

"Abstract type representing bedrock and ice elevation"
abstract type AbstractElevation{T<:Real} end

"Axis-aligned bounding box"
struct AABB{T<:Real}
    xmin::T; xmax::T
    ymin::T; ymax::T
    zmin::T; zmax::T
end

"Construct AABB from coordinates"
AABB(xs,ys,zs) = AABB(extrema(xs)...,extrema(ys)...,extrema(zs)...)

"Elevation data on grid"
struct DataElevation{T, M<:AbstractMatrix{T}} <: AbstractElevation{T}
    x::M; y::M; z_bed::M; z_surf::M
    rotation::M
    #domain::AABB{T}
    #rotated_domain::AABB{T}
end

"Load elevation data from HDF5 file."
function load_elevation(path::AbstractString)
    fid    = h5open(path, "r")
    x      = read(fid,"glacier/x")
    y      = read(fid,"glacier/y")
    z_bed  = read(fid,"glacier/z_bed")
    z_surf = read(fid,"glacier/z_surf")
    R      = read(fid,"glacier/R")
    close(fid)
    return DataElevation(x,y,z_bed,z_surf,R)
end

"Assert that the sampling of a raster is `Points`, see README.  Returns the raster."
function assert_Point(ra::Union{AbstractRaster, AbstractRasterStack})
    out = all(isa.(DimensionalData.Dimensions.LookupArrays.sampling(ra),
                     DimensionalData.Dimensions.LookupArrays.Points))
    @assert out "Raster has not Point sampling!"
    return ra
end
