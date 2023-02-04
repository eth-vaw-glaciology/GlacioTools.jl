using GlacioTools, Test

# Downloading
##############

@testset begin "Download from folder"
    tmpdir = mktempdir() # cleaned up automatically
    GlacioTools.download_file(*("file://", @__DIR__, "/testfiles/data_PlaineMorte.h5"), tmpdir)
    @test isfile(joinpath(tmpdir, "data_PlaineMorte.h5"))

    tmpdir = mktempdir() # cleaned up automatically
    download_dir = joinpath(@__DIR__,"testfiles/")
    get_all_data(download_dir, tmpdir)
    @test isfile(joinpath(tmpdir, "data_PlaineMorte.h5")) && isfile(joinpath(tmpdir, "test.jld2"))
    @test isdir(joinpath(tmpdir, "zip-folder"))
    @test isfile(joinpath(tmpdir, "zip-folder/data_PlaineMorte.h5")) && isfile(joinpath(tmpdir, "zip-folder/test.jld2"))
    # test forcing download
    t1 = stat(joinpath(tmpdir, "data_PlaineMorte.h5")).mtime
    get_all_data(download_dir, tmpdir; preproc=false)
    @test t1 == stat(joinpath(tmpdir, "data_PlaineMorte.h5")).mtime
    get_all_data(download_dir, tmpdir; preproc=false, force_download=true)
    @test t1 < stat(joinpath(tmpdir, "data_PlaineMorte.h5")).mtime
end
@testset begin "Download from URL"
    # download from URL
    tmpdir = mktempdir() # cleaned up automatically
    GlacioTools.download_file("https://github.com/pohlan/GlacioTools.jl/raw/main/test/testfiles/test.jld2", tmpdir)
    @test isfile(joinpath(tmpdir, "test.jld2"))

    tmpdir = mktempdir() # cleaned up automatically
    get_all_data("https://github.com/pohlan/GlacioTools.jl/raw/main/test/testfiles/test.jld2", tmpdir)
    @test isfile(joinpath(tmpdir, "test.jld2"))
end
@testset begin "Download with Dict"
    tmpdir = mktempdir() # cleaned up automatically
    di = Dict(:a => "https://github.com/pohlan/GlacioTools.jl/raw/main/test/testfiles/test.jld2",
              :b => ["https://github.com/pohlan/GlacioTools.jl/raw/main/test/testfiles/test.zip", "https://github.com/pohlan/GlacioTools.jl/raw/main/test/testfiles/data_PlaineMorte.h5"]
              )
    get_all_data(di, tmpdir)
    @test isfile(joinpath(tmpdir, "data_PlaineMorte.h5")) && isfile(joinpath(tmpdir, "test.jld2"))
    @test isdir(joinpath(tmpdir, "zip-folder"))
    @test isfile(joinpath(tmpdir, "zip-folder/data_PlaineMorte.h5")) && isfile(joinpath(tmpdir, "zip-folder/test.jld2"))
end
@testset begin "signed_distance"
    Point2 = GlacioTools.Point2
    p = Point2(1,1)
    poly = [Point2(0,0), Point2(0,2), Point2(2,2), Point2(2,0)]
    @test GlacioTools.signed_distance(p, poly)==-1
    p = Point2(1,0.5)
    @test GlacioTools.signed_distance(p, poly)==-0.5
end
