using GlacioTools, Test

# Downloading
##############

@testset begin "Download from folder"
    tmpdir = mktempdir() # cleaned up automatically
    download_dir = joinpath(@__DIR__,"testfiles/")
    # download from folder
    get_all_data(download_dir, tmpdir)
    @test isfile(joinpath(tmpdir, "data_PlaineMorte.h5")) && isfile(joinpath(tmpdir, "test.jld2"))
    @test 1e-3 > @elapsed get_all_data(download_dir, tmpdir)                     # test force_download = false
    @test 1e-3 < @elapsed get_all_data(download_dir, tmpdir; force_download=true) # test force_download = true
end
@testset begin "Download from URL"
    # download from URL
    tmpdir = mktempdir() # cleaned up automatically
    get_all_data("https://github.com/pohlan/GlacioTools.jl/raw/main/test/testfiles/test.jld2", tmpdir)
    @test isfile(joinpath(tmpdir, "test.jld2"))

    # test the fetch_glacier routine without download
    # a = fetch_glacier("PlaineMorte","A55f/03",datadir="test/testfiles/")
    # @test size(a.x) == (531,257)
end
