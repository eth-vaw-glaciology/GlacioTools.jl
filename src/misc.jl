
"""
    get_plot_extent(p)


Get extent of a Plots-plot as tuple ((x1,x2), (y1,y2)).

Ref: https://github.com/JuliaPlots/Plots.jl/issues/394
"""
get_plot_extent(p) = round.(Int, (p[1].o.get_xlim())), round.(Int, (p[1].o.get_ylim()))

"""
    get_plot_center(p)


Get center of a Plots-plot.

Ref: https://github.com/JuliaPlots/Plots.jl/issues/394
"""
get_plot_center(p) = round.(Int, (mean(p[1].o.get_xlim()), mean(p[1].o.get_ylim())) )


"""
    boxcar(A::AbstractArray, window::AbstractArray, weights)

Moving average filter with spatially dependent filter window and weighting of cells.
Filters over ±window.
"""
function boxcar(A::AbstractArray{T,N}, window::AbstractArray{<:Integer,N},
                weights::AbstractArray{<:Number,N}=ones(size(A)...)) where {T,N}
    out = similar(A)
    R = CartesianIndices(size(A))
    I1, Iend = first(R), last(R)
    Threads.@threads for I in R # @inbounds does not help
        out[I] = NaN
        n, s = 0, zero(eltype(out))
        I_ul = CartesianIndex(I1.I.*window[I])
        for J in CartesianIndices(UnitRange.(max(I1, I-I_ul).I , min(Iend, I+I_ul).I) )
            # used to be CartesianRange(max(I1, I-I_l), min(Iend, I+I_u) )
            # now it is probably something simpler than what I use above
            AJ, w = A[J], weights[J]
            if !isnan(AJ) && w!=0
                s += A[J]
                n += 1
            end
        end
        if n==0
            out[I] = A[I]
        else
            out[I] = s/n # note: ==NaN if n==s==0
        end
    end
    out
end


"""
    smooth_surface(x, y, surfdem, beddem, icethicknesses, mask= surfdem.>=beddem; minwindow=0)

Smooth the surface DEM of ±`icethicknesses`.  Only smoothes where
there is ice and only using ice-covered cells.
"""
function smooth_surface(x, y, surfdem, beddem, icethicknesses, mask=surfdem.>=beddem; minwindow=0)
    dx = x[2]-x[1]
    @assert y[2]-y[1]==dx
    thickness = surfdem .- beddem
    thickness[isnan.(thickness)] .= 0
    window = round.(Int, icethicknesses .* thickness ./ dx)
    window[window.<minwindow] .= minwindow
    surface_smooth = boxcar(surfdem, window, mask)
    return surface_smooth
end


"""
    smooth_surface(ra::Raster(Stack), icethicknesses; minwindow=0)

Smooth the surface DEM of ±`icethicknesses`.  Only smoothes where
there is ice and only using ice-covered cells.
"""
function smooth_surface(ra::RasterStack, icethicknesses; minwindow=0)
    bed, surface = ra[:bed], ra[:surface]
    x,y = get_x(bed), get_y(bed)
    dx = get_dx(bed)
    thickness = surface .- bed
    return smooth_surface(x, y, surface, bed, icethicknesses, ra[:rmask]; minwindow)
end
function smooth_surface_add_stack(ra::RasterStack, icethicknesses; minwindow=0)
    ss = smooth_surface(ra, icethicknesses; minwindow=minwindow)
    # recreate geostack
    return Rasters.DimensionalData.rebuild_from_arrays(ra, (;NamedTuple(ra)..., surface=ss, surface_rough=ra[:surface]))
end
