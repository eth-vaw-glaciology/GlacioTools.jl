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
   # @assert y[2]-y[1]==dx
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

using GeometryBasics: Point2

"""
    signed_distance(p::Point2{T}, poly::AbstractVector{Point2{T}}) where {T}

Returns the distance of a point `p` to the polygon `poly`.
"""
function signed_distance(p::Point2, poly::AbstractVector{<:Point2})
    if poly[1]==poly[end]
        poly = poly[1:end-1]
    end
    d = dot(p - poly[1], p - poly[1])
    s = 1.0
    j = length(poly)
    for i in eachindex(poly)
        e = poly[j] - poly[i]
        w = p - poly[i]
        b = w - e .* clamp(dot(w, e) / dot(e, e), 0.0, 1.0)
        d = min(d, dot(b, b))
        c = p[2] >= poly[i][2], p[2] < poly[j][2], e[1] * w[2] > e[2] * w[1]
        if all(c) || all(.!c)
            s = -s
        end
        j = i
    end
    return s * sqrt(d)
end

# from WhereTheWaterFlows
"Return CartesianIndices corresponding to the 8 neighbors and the point itself."
function iterate_D9(I::CartesianIndex, ar::AbstractMatrix)
    R = CartesianIndices(ar)
    one = CartesianIndex(1,1)
    I1, Iend = first(R), last(R)
    return max(I1, I-one):min(Iend, I+one)
end

# from WhereTheWaterFlows
"Return CartesianIndices corresponding to the 24 neighbors and the point itself."
function iterate_D25(I::CartesianIndex, ar::AbstractMatrix)
    R = CartesianIndices(ar)
    two = CartesianIndex(2,2)
    I1, Iend = first(R), last(R)
    return max(I1, I-two):min(Iend, I+two)
end


"""
    iterate_D5(I::CartesianIndex, ar::AbstractMatrix)

Return CartesianIndices corresponding to the 4 non-diagional neighbors and the point itself.

TODO: this is pretty slow, about 5x slower than iterate_D9
"""
function iterate_D5(I::CartesianIndex, ar::AbstractMatrix)
    i, j = Tuple(I)
    R = CartesianIndices(ar)
    I1, Iend = first(R), last(R)

    out = (CartesianIndex(i,j), )
    for it in (CartesianIndex(i-1,j),
               CartesianIndex(i+1,j),
               CartesianIndex(i,j-1),
               CartesianIndex(i,j+1))
        if it in max(I1, I-I1):min(Iend, I+I1)
            out = (it, out...)
        end
    end
    return out
end


"""
    mask_contiguous(mask, IJ, out=similar(mask, Bool).*false )

Create a new mask (or update the optional 3rd argument) which
- marks all points connected to IJ and
- where each masked point has the same value as mask[IJ]
"""
function mask_contiguous(mask, IJ, out=similar(mask, Bool).*false )
    val = mask[IJ]
    queue = [IJ]
    out[IJ] = true
    while length(queue)>0
        # get first element
        IJ = popfirst!(queue)
        # push onto queue
        for ij in iterate_D9(IJ, out)
            out[ij] && continue # already done or queued
            if mask[ij]==val
                out[ij] = true
                push!(queue, ij) # queue new point
            end
        end
        # make sure it does not blow up
        length(queue)>prod(size(mask)) && error("oh no, queue got too big")
    end
    return out
end

"""
    get_boudary_cells(mask, value)

Get all the cells of mask of value `value` which border
on cells of different value.
"""
function get_boudary_cells(mask, value)
    out = similar(mask, Bool) .* false

    for IJ in CartesianIndices(mask)
        mask[IJ]!=value && continue # don't process
        for ij in iterate_D9(IJ, mask)
            ij==IJ && continue # don't process the point itself
            if mask[ij]!=value # found the boundary
                out[IJ] = true
                continue
            end
        end
    end
    return out
end

"""
     dilate(img::AbstractArray{<:Bool})

Image processing "dilation", in short, makes regions of "true" larger by a pixel;
opposite of "erosion".
https://juliaimages.org/stable/examples/image_morphology/image_morphology/#Dilation
"""
function dilate(img::AbstractArray{<:Bool})
    out = img .* false
    for IJ in CartesianIndices(img)
        # use a D9 Kernel as is standard
        for ij in iterate_D9(IJ, img)
            if img[ij] # found connection
                out[IJ] = true
                break
            end
        end
    end
    return out
end


"""
    piecewiselinear(xs::AbstractVector, ys::AbstractVector)

Return a piecewise linear function `f(x) -> y`` through points (xs,ys).
Extrapolation with x<xs[1]->ys[1] and for x>xs[end] -> ys[end] .

Notes:
- xs needs monotonic, i.e. to be sorted in ascending or decending order
- when length(xs)==1, then it returns a constant function.
- currently extrapolates using the last value. (TODO)

Almost 2x faster than using Interpolations.jl.
"""
function piecewiselinear(xs::AbstractVector, ys::AbstractVector)
    if length(xs)==1
        @assert length(xs)==length(ys)
        return xq -> ys[1]
    end
    rats = diff(ys)./diff(xs)
    if issorted(xs)
        return function (xq)
            # xq<xs[1] && error("cannot extrapolate")
            # xq>xs[end] && error("cannot extrapolate")
            # xq==xs[end] && return ys[end]
            xq<=xs[1] && return ys[1]
            xq>=xs[end] && return ys[end]
            ii = searchsortedlast(xs, xq)
            out = ys[ii] - (xs[ii]-xq)*rats[ii]
            return out
        end
    elseif issorted(xs, Base.Order.Reverse)
        return function (xq)
            # xq>xs[1] && error("cannot extrapolate")
            # xq<xs[end] && error("cannot extrapolate")
            # xq==xs[end] && return ys[end]
            xq>=xs[1] && return ys[1]
            xq<=xs[end] && return ys[end]
            ii = searchsortedlast(xs, xq, Base.Order.Reverse)
            out = ys[ii] - (xs[ii]-xq)*rats[ii]
            return out
        end
    else
        error("xs must be sorted")
    end
end

"""
    parameterized_curve(x,y)

Return curve-length parameterized curve-function and
length range.

Note:
- only works for monotonic `x`!

TODO: make it work for non-monotoic `x``.

Example
```
julia> x,y = [1,2,3], [6, 11, 2];

julia> pf, span = parameterized_curve(x,y);

julia> xy = hcat(pf.(r)...); plot(xy[1,:], xy[2,:])
````
"""
function parameterized_curve(x,y)
    # TODO: make dimension agnostic sometime
    pathlength = zeros(typeof(x[1]), length(x)-1) # assume all same datatype
    for c in (x,y)
        pathlength .+= diff(c).^2
    end
    pathlength = [0, cumsum(sqrt.(pathlength))...]
    fx = piecewiselinear(pathlength, x)
    fy = piecewiselinear(pathlength, y)
    return pl -> [fx(pl), fy(pl)], (0, pathlength[end])
end
