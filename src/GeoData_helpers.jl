get_dx(ga::Union{GeoArray,GeoStack}) = diff(val(dims(ga)[1])[1:2])[1]
get_dy(ga::Union{GeoArray,GeoStack}) = diff(val(dims(ga)[2])[1:2])[1]
get_x(ga::Union{GeoArray,GeoStack})  = val(dims(ga)[1])
get_y(ga::Union{GeoArray,GeoStack})  = val(dims(ga)[2])


"""
    missing2nan(ar, T=Float32)
Convert `missing` or `==missingval` to NaN and return a `GeoArray{Float32,2}`.
"""
function missing2nan(ar::A, T=Float32) where A<:AbstractGeoArray{<:Number,2}
    # convert to T first
    data = convert(Matrix{T}, ar.data)
    data[data.==ar.missingval] .= NaN
    return GeoArray(data; ar.dims, ar.name, ar.refdims, ar.metadata)
end
function missing2nan(ar::AbstractGeoArray, T=Float32) # Union{Missing,...}
    data = convert(Matrix{Union{T,Missing}}, copy(ar.data))
    data[ismissing.(data)] .= NaN
    data = convert(Matrix{T}, data)
    return GeoArray(data; ar.dims, ar.name, ar.refdims, ar.metadata)
end
