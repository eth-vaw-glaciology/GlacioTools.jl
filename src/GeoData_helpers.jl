# From FourDAntarcticaSubglacialRouting.jl #
############################################
get_dx(ga::Union{Raster,RasterStack}) = diff(val(dims(ga)[1])[1:2])[1]
get_dy(ga::Union{Raster,RasterStack}) = diff(val(dims(ga)[2])[1:2])[1]
get_x(ga::Union{Raster,RasterStack})  = val(dims(ga)[1])
get_y(ga::Union{Raster,RasterStack})  = val(dims(ga)[2])

"""
    Box(x::Tuple, y::Tuple)
Get an array that can be used to crop a Raster to specified x and y coordinate ranges.
"""
Box(x::Tuple, y::Tuple) = (X(x[1]..x[2]), Y(y[1]..y[2]))

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
