using GlacioTools, Test
const GT = GlacioTools

# Not included because it downloads tons of data:
#include("runtest_antarctica.jl")

include("runtest_misc_filereaders.jl")


# Downloading
##############

@testset begin "Download from folder"
    tmpdir = mktempdir() # cleaned up automatically
    GT.download_file(*("file://", @__DIR__, "/testfiles/data_PlaineMorte.h5"), tmpdir)
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
    GT.download_file("https://github.com/pohlan/GlacioTools.jl/raw/main/test/testfiles/test.jld2", tmpdir)
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
    Point2 = GT.Point2
    p = Point2(1,1)
    poly = [Point2(0,0), Point2(0,2), Point2(2,2), Point2(2,0)]
    @test GT.signed_distance(p, poly)==-1
    p = Point2(1,0.5)
    @test GT.signed_distance(p, poly)==-0.5

    p = Point2(-2.355e6, 1.247e6)
    poly = Point2{Float64}[[-2.3959741684207064e6, 1.2260411126201542e6], [-2.3921964087091144e6, 1.2279839604718303e6], [-2.3871538317031735e6, 1.2316512892034235e6], [-2.3834934311744734e6, 1.2332374627658604e6], [-2.3825696707886816e6, 1.2336377589330368e6], [-2.3789023420570884e6, 1.2350130072073843e6], [-2.3773742884189244e6, 1.2389859466666104e6], [-2.3753878186893114e6, 1.2405140003047744e6], [-2.3718732953215344e6, 1.2414308324876728e6], [-2.3683587719537574e6, 1.241889248579122e6], [-2.3659138861326952e6, 1.2437229129449185e6], [-2.3645717426209664e6, 1.2451489404261303e6], [-2.3642162076957542e6, 1.2445464715182958e6], [-2.3639783297631517e6, 1.2437196253065444e6], [-2.3640980837595845e6, 1.2425433346426988e6], [-2.364553105408363e6, 1.2415822513549794e6], [-2.365728046013621e6, 1.2391666808686983e6], [-2.365462220282019e6, 1.237836131785414e6], [-2.36528887242531e6, 1.2365651168825321e6], [-2.3648623203496826e6, 1.2363385081741505e6], [-2.3640074137365837e6, 1.2364080732286812e6], [-2.3631152107814476e6, 1.2367223928205357e6], [-2.3625412163186297e6, 1.2368297928570386e6], [-2.3619569852590933e6, 1.2365051047234589e6], [-2.3613836239213212e6, 1.2363071030245996e6], [-2.360558641514598e6, 1.2361903821204058e6], [-2.3601660447904244e6, 1.2360110154107208e6], [-2.35978034191459e6, 1.2356002219629642e6], [-2.359579921387345e6, 1.234799975375218e6], [-2.358992716131619e6, 1.2342141972526505e6], [-2.3587923642756757e6, 1.233960551989749e6], [-2.3588210724558877e6, 1.232970010774121e6], [-2.3586244135574363e6, 1.2325055468345066e6], [-2.3577965573904524e6, 1.2322791957180873e6], [-2.3577116866371357e6, 1.2321183526413098e6], [-2.359095619319546e6, 1.2290121926207894e6], [-2.3607615489188596e6, 1.2252836835175639e6], [-2.362665468460932e6, 1.2233400989850312e6], [-2.3638554181747274e6, 1.2224278042044549e6], [-2.3657990027072597e6, 1.221317184471579e6], [-2.367068282401975e6, 1.2211188595192798e6], [-2.372700711047273e6, 1.2179456602824922e6], [-2.3748426205321047e6, 1.2161210707213392e6], [-2.3772225199596956e6, 1.2129478714845516e6], [-2.380395719196483e6, 1.2073154428392532e6], [-2.3798404093300453e6, 1.2058081732017791e6], [-2.3784917996544107e6, 1.2029522938886702e6], [-2.378412469673491e6, 1.200651724441999e6], [-2.3780158197688926e6, 1.199382444747284e6], [-2.37865045961625e6, 1.1964472354532555e6], [-2.378412469673491e6, 1.192480736407271e6], [-2.379681749368206e6, 1.1909734667697968e6], [-2.382458298700395e6, 1.189624857094162e6], [-2.386583457708219e6, 1.1858170180100168e6], [-2.3876940774410944e6, 1.1833577886015063e6], [-2.3883978600860313e6, 1.182115819228089e6], [-2.3884873772502933e6, 1.1829214737064478e6], [-2.389597996983169e6, 1.1870466327142718e6], [-2.390311966811446e6, 1.1884745723708263e6], [-2.3930885161436354e6, 1.1924410714168109e6], [-2.395389085590306e6, 1.1952969507299198e6], [-2.3972136751514594e6, 1.1964075704627954e6], [-2.398562284827094e6, 1.1981528300430288e6], [-2.39967290455997e6, 1.1992634497759044e6], [-2.4015768241020422e6, 1.1997394296614225e6], [-2.403480743644115e6, 1.2010087093561376e6], [-2.4043533734342316e6, 1.2026746389554513e6], [-2.4054639931671075e6, 1.2036265987264875e6], [-2.4068919328236617e6, 1.2035472687455679e6], [-2.4077645626137783e6, 1.203785258688327e6], [-2.4087165223848145e6, 1.2052131983448814e6], [-2.409351162232172e6, 1.2067204679823555e6], [-2.410461781965048e6, 1.207910417696151e6], [-2.411572401697924e6, 1.2083863975816693e6], [-2.4095098221940096e6, 1.208426062572129e6], [-2.4075265726710176e6, 1.2082674026102896e6], [-2.403560073625033e6, 1.208664052514888e6], [-2.3990382647126103e6, 1.2099333322096032e6], [-2.3965790353041e6, 1.2144551411220257e6], [-2.394199135876509e6, 1.2197702498436451e6], [-2.3953890855903043e6, 1.2242127287751478e6], [-2.3959741684207064e6, 1.2260411126201542e6]]
    @test !isnan(GT.signed_distance(p, poly))
end


@testset begin "mask_contiguous"
    mask = ones(Int, 30, 31)
    mask[[CartesianIndex(i,j) for i=1:10, j=1:15 if i<j]] .= 3

    @test sum(GT.mask_contiguous(mask, CartesianIndex(4,5))) == sum(mask.==3)
    @test sum(GT.mask_contiguous(mask, CartesianIndex(5,2))) == sum(mask.==1)

    mask2 = copy(mask)
    mask2[[CartesianIndex(i,j) for i=20:27, j=17:20 if i<j]] .= 3

    @test sum(GT.mask_contiguous(mask2, CartesianIndex(4,5))) == sum(mask.==3)
    @test sum(GT.mask_contiguous(mask2, CartesianIndex(5,2))) == sum(mask.==1)
end
