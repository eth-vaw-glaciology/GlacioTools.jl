# From FourDAntarcticaSubglacialRouting.jl #
############################################
get_x(ra::Union{Raster,RasterStack})  = dims(ra)[1].val
get_y(ra::Union{Raster,RasterStack})  = dims(ra)[2].val
get_dx(ra::Union{Raster,RasterStack}) = diff(get_x(ra)[1:2])[1]
get_dy(ra::Union{Raster,RasterStack}) = diff(get_y(ra)[1:2][1:2])[1]

"""
    Box(x::Tuple, y::Tuple)

Get an index-tuple that can be used to crop a Raster to specified x and y coordinate ranges.

It is really just this:

    (X(x[1]..x[2]), Y(y[1]..y[2]))
"""
Box(x::Tuple, y::Tuple) = (X(x[1]..x[2]), Y(y[1]..y[2]))

## Geography
function inbox(pt, box)
    xr,yr = val.(val.(box))
    (xr[1]<=pt[1]<=xr[2]) && (yr[1]<=pt[2]<=yr[2])
end
center(box) = sum.(val.(val.(box))).÷2

"""
    box2cartesian(box, ra::Raster)

Translates a coordinate box into CartesianIndices.
"""
box2cartesian(box, ra::Raster) = CartesianIndices(Rasters.DimensionalData.Dimensions.dims2indices(ra, box))


# From FastIce.jl/GeoData #
###########################

"Helper function to mask, trim and pad bedrock and ice thickness data given a glacier polygon."
@views mask_trim(rasterDat, poly, pad) = trim(mask(rasterDat; with=poly); pad=pad)

"Filter out all values of `A` based on `mask`."
function my_filter(A, mask)
    return [A[i] for i in eachindex(A) if mask[i] != 0]
end

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
