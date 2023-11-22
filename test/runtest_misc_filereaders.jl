@testset "ASCII-Grid bin files" begin
    da = GlacioTools.read_agr("testfiles/ascii-grid.bin")
    @test size(da) == (527, 369)
    @test isequal(sum(da), NaN32)
    @test sum(da[.!isnan.(da)]) ≈ 7.549487f7
end
