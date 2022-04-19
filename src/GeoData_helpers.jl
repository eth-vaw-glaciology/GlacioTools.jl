get_dx(ga::Union{Raster,RasterStack}) = diff(val(dims(ga)[1])[1:2])[1]
get_dy(ga::Union{Raster,RasterStack}) = diff(val(dims(ga)[2])[1:2])[1]
get_x(ga::Union{Raster,RasterStack})  = val(dims(ga)[1])
get_y(ga::Union{Raster,RasterStack})  = val(dims(ga)[2])

"""
    Box(x::Tuple, y::Tuple)
Get an array that can be used to crop a Raster to specified x and y coordinate ranges.
"""
Box(x::Tuple, y::Tuple) = (X(x[1]..x[2]), Y(y[1]..y[2]))

"""
    missing2nan(ar, T=Float32)
Convert `missing` or `==missingval` to NaN and return a `Raster{Float32,2}`.
"""
function missing2nan(ar::A, T=Float32) where A<:AbstractRaster{<:Number,2}
    # convert to T first
    data = convert(Matrix{T}, ar.data)
    data[data.==ar.missingval] .= NaN
    return Raster(data; ar.dims, ar.name, ar.refdims, ar.metadata)
end
function missing2nan(ar::AbstractRaster, T=Float32) # Union{Missing,...}
    data = convert(Matrix{Union{T,Missing}}, copy(ar.data))
    data[ismissing.(data)] .= NaN
    data = convert(Matrix{T}, data)
    return Raster(data; ar.dims, ar.name, ar.refdims, ar.metadata)
end
