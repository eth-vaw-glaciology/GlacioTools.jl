using GlacioTools, Test

@testset begin
    tmp = "test/tempdir/"
    # download from folder
    get_all_data("test/testfiles",tmp)
    @test isfile("test/tempdir/alps/data_PlaineMorte.h5") && isfile("test/tempdir/test.jld2")
    @test 1e-3 > @elapsed get_all_data("test/testfiles",tmp)                     # test force_download = false
    @test 1e-3 < @elapsed get_all_data("test/testfiles",tmp;force_download=true) # test force_download = true
    run(`rm -r $tmp`)

    # download from URL
    get_all_data("https://github.com/pohlan/GlacioTools.jl/raw/main/test/testfiles/test.jld2",tmp)
    @test isfile(tmp*"test.jld2")
    run(`rm -r $tmp`)

    # test the fetch_glacier routine without download
    # a = fetch_glacier("PlaineMorte","A55f/03",datadir="test/testfiles/")
    # @test size(a.x) == (531,257)
end
