# Readers of special (some in-house) file types
#
# TODO:
# - ascci writer
# - xyz reader / writer

"""
Holds a rectangular grid, mirrors the ASCII-grid .agr, .grid, .asc, .bin files.

Note that the indexing into the value matrix ins awkward, see Ref
below.

Ref: https://en.wikipedia.org/wiki/Esri_grid
"""
struct AGR{T} # Ascii GRid
    v::Matrix{T}    # values: orientation is awkward, as in the AGR file, see https://en.wikipedia.org/wiki/Esri_grid
    nc::Int64       # NCOLS TODO: remove those
    nr::Int64       # NROWS
    xll::Float64    # XLLCORNER
    yll::Float64    # YLLCORNER
    dx::Float64     # CELLSIZE
    NA::T           # NODATA
    hasutm::Bool    # Matthias sometimes abuses the NODATA_value field as UTM-zone field
    extra_header::Vector{Float32} # the .bin files have space for
                                  # extra information in the header
    function AGR{T}(va, nc, nr, xll, yll, dx, NA,
                    hasutm=false, extra_header=zeros(Float32,6)) where T
        @assert size(va)==(nr,nc)
        @assert dx>=0
        new{T}(va, nc, nr, xll, yll, dx, NA, hasutm, extra_header)
    end
end
AGR(va::Matrix{T}, nc, nr, xll, yll, dx, NA, hasutm=false, extra_header=zeros(Float32,6)) where {T} = AGR{T}(va, nc, nr, xll, yll, dx, NA, hasutm, extra_header)


Base.size(g::AGR) = size(g.v)
Base.length(g::AGR) = length(g.v)
function Base.:(==)(g1::AGR,g2::AGR)
    g1.v==g2.v &&
    g1.nc==g2.nc &&
    g1.nr==g2.nr &&
    g1.xll==g2.xll &&
    g1.yll==g2.yll &&
    g1.dx==g2.dx &&
    g1.NA==g2.NA &&
    g1.extra_header==g2.extra_header
end

"""
    agr2raster(agr::AGR{T}; NA=convert(T,NaN)) where T

Transform asciigrid to raster.
"""
function agr2raster(agr::AGR{T}; NA=convert(T,NaN)) where T
    #    v = my_rotr90(agr.v)
    v = rotr90(agr.v)
    if agr.hasutm
        proj = "+proj=utm +zone=$(Int(agr.NA)) +datum=WSG84"
    else
        # only swap NA if it has a FILL value:
        if !isequal(NA, agr.NA)
            for i=eachindex(v)
                if isequal(v[i], agr.NA); v[i] = NA end
            end
        end
        proj = ""
    end
    return Raster(v, (X(range(agr.xll+agr.dx/2, step=agr.dx, length=agr.nc)),
                      Y(range(agr.yll+agr.dx/2, step=agr.dx, length=agr.nr))),
                     crs=proj)
end

# AGR file readers
##################

"""
Check if .bin file nor not.  With IO types this may not work 100%.
"""
isbin_file(fn::AbstractString) = splitext(fn)[2]==".bin"
function isbin_file(fn::IO)
    if :name in fieldnames(typeof(fn))
        f = fn.name
    else # this should catch TranscodingStreams.jl
        f = fn.stream.name
    end
    f = strip(f, ['<','>']) # often looks like "<file asdf.ext>"
    f = strip(f)
    isbin_file(f)
end

"""
    read_agr(fl, T=Float32; NA=convert(T,NaN))

Read a Ascii grid https://en.wikipedia.org/wiki/Esri_grid.  Can also
read from .gz compressed files, if `import CodecZlib`.

In:
- fn -- file name

Optional keywords:
- T -- type of output array. Defaults to Float32.
- NA -- replace the fill value with this. Defaults to use NaN.

"""
function read_agr(fl, T=Float32; NA=convert(T,NaN))
    agr = _read_agr(fl, T)
    (;v, hasutm) = agr
    NA_old = agr.NA
    # replace missing value with something else
    if !hasutm && NA_old!=NA
        _refill!(v, NA_old, NA)
    end
    return agr2raster(agr; NA)
end

function _read_agr(fl::AbstractString, T=Float32)
    if !isfile(fl)
        error("File $fl does not exist!")
    end
    if endswith(fl, ".gz")
        @eval import CodecZlib
        io = CodecZlib.GzipDecompressorStream(open(fl))
    else
        io = open(fl, "r")
    end
    out = _read_agr(io, T)
    close(io)
    return out
end

function _read_agr(io::IO, T=Float32)
    if isbin_file(io)
        toT = (io,T) -> convert(T, read(io, Float32))
    else
        toT = (io,T) -> parse(T, split(readline(io))[2])
    end

    local va::Matrix{T}
    # read header
    nc = toT(io, Int)
    nr = toT(io, Int)
    xll = toT(io,Float64)
    yll = toT(io,Float64)
    dx = toT(io,Float64)
    # Matthias sometimes abuses the NODATA_value field as UTM-zone field
    # in non-binary grids:
    if isbin_file(io)
        fill = toT(io,T)
        hasutm = false
    else
        prop, val = split(readline(io))
        if lowercase(prop)=="nodata_value"
            hasutm = false
        elseif lowercase(prop)=="utm_zone"
            hasutm = true
        else
            error("Unrecognized last header field: $prop")
        end
        fill = parse(T,val)
    end

    if isbin_file(io)
        # read extra header
        extra_header = Array{Float32}(undef, 6)
        read!(io, extra_header)
        # read values
        va = Array{Float32}(undef, nc, nr)
        va = permutedims(read!(io, va))
        if eltype(va)!=T
            error("Not implemented yet")
        end
        if !eof(io)
            @warn("End-of-file was not reached!")
        end
    else
        va = Array{T}(undef, nr, nc)
        # no extra header for ascii .agr
        extra_header = zeros(Float32, 6)
        # read values
        tmp = split(read(io, String))
        if length(tmp)!=nc*nr
            error("Something's wrong with the file/stream $io.  It should contain $(nc*nr) values but has $(length(tmp))")
        end
        for i=1:nr, j=1:nc
            va[i,j] = parse(T, tmp[(i-1)*nc + j])
        end
    end
    # make a AGR data structure
    AGR(va, nc, nr, xll, yll, dx, fill, hasutm, extra_header)
end
# A function barrier is needed, thus this helper function:
function _refill!(a::Array, oldfill, newfill)
    @inbounds for i in eachindex(a)
        if a[i]==oldfill
            a[i] = newfill
        end
    end
    newfill
end

"""
    get_utm_asciigrid(fl_or_io)

Get the UTM zone of an ASCII-grid file (a unofficial format change used by Matthias).

If there is no zone return "".
"""
function get_utm_asciigrid(fl)
    open(fl) do io
        get_utm_asciigrid(io)
    end
end
function get_utm_asciigrid(io::IO)
    T = Float32
    if isbin_file(io)
        return ""
        # error("Ascii-bin files cannot have UTM")
    else
        toT = (io,T) -> parse(T, split(readline(io))[2])
    end

    # read header
    _ = toT(io, Int)
    _ = toT(io, Int)
    _ = toT(io,T)
    _ = toT(io,T)
    _ = toT(io,T)
    # Matthias sometimes abuses the NODATA_value field as UTM-zone field
    # in non-binary grids:
    prop, val = split(readline(io))
    if lowercase(prop)=="nodata_value"
        return ""
        # error("No field UTM_ZONE found")
    end
    utm = parse(Int,val)
    return "+proj=utm +zone=$(utm) +datum=WGS84"
end

# """
#     write_agr(g::AGR{T_}, fn::AbstractString; NA=nothing, T=T_) where T_

# Write Ascii grid.

# Output format is determined by the extension:
# - .bin -- binary
# - .agr -- ASCII
# - .grid or .asc are renamed to .agr

# Always writes .bin in Float32

# Optional:
# - T: choose output type: Int or Float (only for ASCII files)
# """
# function write_agr(g::AGR{T_}, fn::AbstractString; NA=nothing, T=T_) where T_
#     if !isbin_file(fn)
#         ext = splitext(fn)[2]
#         if ext==".grid" || ext==".asc"
#             # println("Changing extension to .agr")
#             # ext = ".agr"
#         elseif ext!=".agr"
#             error("unsupported extension")
#         end
#     end
#     if isbin_file(fn)
#         convInt(x) = convert(Float32, x)
#         convT(x) = convert(Float32, x)
#         convMT = convT
#     else # conversion to ASCII
#         convInt = (x) -> @sprintf("%d\n", x)
#         if T<:AbstractFloat
#             convT = (x) ->  @sprintf("%f\n", x)
#             convMT = (x) ->  @sprintf("%f  ", x)
#         elseif T<:Integer
#             convT = (x) ->  @sprintf("%d\n", x)
#             convMT = (x) ->  @sprintf("%d  ", x)
#         end
#     end
#     open(fn, "w") do io
#         # write header
#         isbin_file(fn) || write(io, "ncols         ")
#         write(io, convInt(g.nc))
#         isbin_file(fn) || write(io, "nrows         ")
#         write(io, convInt(g.nr))
#         isbin_file(fn) || write(io, "xllcorner     ")
#         write(io, convT(g.xll))
#         isbin_file(fn) || write(io, "yllcorner     ")
#         write(io, convT(g.yll))
#         isbin_file(fn) || write(io, "cellsize      ")
#         write(io, convT(g.dx))
#         if !isbin_file(fn)
#             if g.hasutm
#                 write(io, "UTM_ZONE      ")
#             else
#                 write(io, "NODATA_value  ")
#             end
#         end
#         # replace NA value if necessary
#         if NA==nothing
#             fill = g.NA
#             va = g.v
#         else
#             fill = NA
#             va = deepcopy(g.v)
#             for i=1:length(va)
#                 if va[i]==g.NA
#                     va[i] = fill
#                 end
#             end
#         end
#         if g.hasutm
#             write(io, convInt(fill))
#         else
#             write(io, convT(fill))
#         end
#         # write extra header if .bin:
#         isbin_file(fn) && write(io, g.extra_header)
#         # write values
#         for i=1:g.nr, j=1:g.nc
#             v = va[i,j]
#             write(io, convMT(v))
#         end
#         return nothing
#     end
# end

# """
# Read .xyn or .xyzn files which can contain several, joined polygons.

# x              y           label
# -2031744.122   833011.310  21
# ...

# x              y           z   label
# -2031744.122   833011.310  23.4  21
# ...

# where the label marks the beginning or end of line (21 --start, 22 --
# inside, 23 -- end).

# There can be embedded polygons, with label sequence 23,20,21.  The
# line with 20 is then dropped.

# The polygon must be closed, i.e. first element equals last.

# Set `hasz` if it also contains a z-coordinate, i.e. a .xyzn file.

# Return:

# - a list of [x,y(,z)] arrays, one for each polygon.
# """
# function read_xyn(fn; hasz=false, fix=false)
#     if !isfile(fn)
#         error("File $fn cannot be found.")
#     end
#     x,y,z,l = open(fn, "r") do io
#         x = Float64[]
#         y = Float64[]
#         z = Float64[]
#         l = Int[]
#         for ls in readlines(io)
#             if hasz
#                 x_,y_,z_,l_ = split(ls)
#             else
#                 x_,y_,l_ = split(ls)
#             end
#             push!(x, parse(Float64, x_))
#             push!(y, parse(Float64, y_))
#             if hasz
#                 push!(z, parse(Float64, z_))
#             end
#             push!(l, parse(Int, l_))
#         end
#         return x,y,z,l
#     end
#     # Split up into individual polygons
#     n = length(x)
#     out = Matrix{Float64}[]
#     is = 1 # start index of one poly
#     iend = -1  # end index of one poly
#     i = 1
#     drop20 = false
#     while i<=n
#         if l[i]!=21
#             error("Malformed file $fn: Expected 21 on line $i.")
#         end
#         # start a new polygon
#         i += 1
#         while true  # go through one poly:
#             if l[i]==23 # end of a poly
#                 if i!=n && l[i+1]==20 # drop the 20-line
#                     drop20 = true
#                 end
#                 break
#             elseif l[i]==22
#                 i+=1
#             else
#                 error("shouldn't get here")
#             end
#         end
#         iend = i
#         if length(is:iend)<3
#             error("Polygon needs at least three vertices")
#         end
#         # fill:
#         if hasz
#             push!(out, (hcat(x[is:iend],y[is:iend],z[is:iend]))' )
#         else
#             push!(out, (hcat(x[is:iend],y[is:iend]))' )
#         end
#         if x[is]!=x[iend] || y[is]!=y[iend]
#             if fix
#                 out[end] = hcat(out[end], out[end][:,1])
#             else
#                 error("Malformed file $fn: polygon $is:$iend not closed")
#             end
#         end

#         is = iend+1
#         i+=1
#         if drop20
#             is+=1
#             i+=1
#             drop20 = false
#         end
#     end
#     return out
# end
# function read_xyz(fn)
#     if !isfile(fn)
#         error("File $fn cannot be found.")
#     end
#     x,y,z = open(fn, "r") do io
#         x = Float64[]
#         y = Float64[]
#         z = Float64[]
#         for ls in readlines(io)
#             x_,y_,z_ = split(ls)
#             push!(x, parse(Float64, x_))
#             push!(y, parse(Float64, y_))
#             push!(z, parse(Float64, z_))
#         end
#         return x,y,z
#     end
#     # # turn into grid
#     tmp = diff(x)
#     dx = median(tmp[tmp.>0])
#     tmp = diff(y)
#     dy = median(tmp[tmp.>0])
#     xrange = extrema(x)
#     yrange = extrema(y)
#     xx = xrange[1]:dx:xrange[2]
#     yy = yrange[1]:dy:yrange[2]
#     zz = zeros(length(xx), length(yy)).*NaN
#     for (X,Y,Z) in zip(x,y,z)
#         i,j = findfirst(isequal(X), xx), findfirst(isequal(Y), yy)
#         zz[i,j] = Z
#     end
#     return (xx,yy,zz)
# end

# """
#     concat_poly(mpoly::Vector)

# Concatenates a split poly, as returned from read_xyn, into one. Also
# returns indices where to split apart again.  The concatenated polygon
# is fully connected, with an edge going back to the first point.  This
# allows to use `inpoly`, at least if the inner polygons have different
# orientation to the outer. Also note that the input and output polygons
# are closed, i.e. last point == first point (this can be fixed by setting
# the option `close_poly=true`)


# Return:
# - bigpoly -- with size==(2,n)
# - splits -- ith poly has indices splits[i]:splits[i+1]-1
# """
# function concat_poly(mpoly::Vector; close_poly=false)
#     T = eltype(mpoly[1])
#     # check that they are all closed
#     if close_poly
#         mpoly = deepcopy(mpoly)
#     end
#     for i=1:length(mpoly)
#         if mpoly[i][:,1]!=mpoly[i][:,end]
#             if close_poly
#                 mpoly[i] = hcat(mpoly[i], mpoly[i][:,1])
#             else
#                 error("All input polys need to be closed.")
#             end
#         end
#     end
#     # total size is sum of sizes plus one extra point for all but the
#     # first poly.
#     totsize = mapreduce(x->size(x,2), +, mpoly) + length(mpoly) -1
#     bigpoly = Array{T}(undef, size(mpoly[1],1), totsize)
#     splits = Int[]
#     is = 1
#     for i=1:length(mpoly)
#         push!(splits, is)
#         bigpoly[:,is:is+size(mpoly[i],2)-1] = mpoly[i]
#         if i==1
#             is = is+size(mpoly[i],2)
#         else
#             # add extra point to connect back to mpoly[1][:,1]
#             is = is+size(mpoly[i],2) + 1
#             bigpoly[:,is-1] = mpoly[1][:,1]
#         end
#     end
#     push!(splits, totsize+1)
#     return bigpoly, splits
# end

# "Finds splitting points in a big-poly"
# function find_poly_splits(bigpoly::Matrix)
#     splits = Int[]
#     p1 = bigpoly[:,1]
#     for i=1:size(bigpoly,2)
#         if p1==bigpoly[:,i]
#             push!(splits,i+1)
#         end
#     end
#     splits
# end

# "Split up concatenated polygon."
# function split_poly(bigpoly::Matrix{T}, splits) where T
#     out = Matrix{T}[]
#     if length(splits)==0
#         return out
#     else
#         push!(out, bigpoly[:,splits[1]:splits[1+1]-1])
#     end
#     for i=2:length(splits)-1
#         # remove the extra point joining to the first poly
#         push!(out, bigpoly[:,splits[i]:splits[i+1]-2])
#     end
#     out
# end
