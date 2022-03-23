"""
    download_file(url, destination_dir;
                    force_download=false,
                    filename=nothing)
Download a file, if it has not been downloaded already.  Uses wget, which needs to64
be installed. Doesn't work for files protected with a password.

Input:
- url -- url for download
- destination_dir -- path of the directory to store the download
Kwargs:
- force_download -- force the download, even if file is present
- filename -- rename file
Output:
- destination_file -- file name (including relative path)
"""
function download_file(url::String, destination_dir::String;
                        force_download=false,
                        filename=nothing)
    mkpath(destination_dir)
    destination_file = if filename === nothing
        joinpath(destination_dir, basename(url))
    else
        joinpath(destination_dir, filename)
    end
    if isfile(destination_file) && !force_download
        # do nothing
        return destination_file
    elseif isfile(destination_file)
        rm(destination_file)
    end
    if startswith(url, "file://")
        # just copy, without the "file://" part
        cp(url[8:end], destination_file, force=force_download)
    else
        Downloads.download(url, destination_file)
    end
    return destination_file
end

"""
    preproc_data(fl, destination_dir)
Unpack downloaded .zip, .gz or .tar files
"""
function preproc_data(fl, destination_dir)
    if splitext(fl)[2]==".zip"
        cd(destination_dir) do
            run(`unzip -ou $fl`)
        end
    elseif splitext(fl)[2]==".gz"
        @assert splitext(splitext(fl)[1])[2]==".tar"
        cd(destination_dir) do
            run(`tar xzf $fl`)
        end
    elseif splitext(fl)[2]==".tar"
        cd(destination_dir) do
            run(`tar xf $fl`)
        end
    end
end

"""
    get_all_data(datas, destionation_dir::String;
                    force_download=false)
Downloads all files in a directory or dictionary.

Input:
- datas -- collection of files that need to be downloaded; can be a dictionary, a folder or just the path of a single file.
- destionation_dir -- path of the directory to store the download
Kwargs:
- force_download -- force the download, even if file is present
"""
function get_all_data(datas::Dict, destionation_dir::String;
                        force_download=false)
    for (k,d) in datas
        print("Downloading $k... ")
        if d isa AbstractString
            try
                fl = download_file(d, destionation_dir; force_download)
                preproc_data(fl, destionation_dir)
            catch e
                println(" ... error: $e")
            end
        else
            for dd in d
                try
                    fl = download_file(dd, destionation_dir; force_download)
                    preproc_data(fl, destionation_dir)
                catch e
                    println(" ... error: $e")
                end
            end
        end
        println("done.")
    end
    nothing
end
function get_all_data(datas::String, destionation_dir::String;
                        force_download=false)
    if isdir(datas)
        d = readdir(datas, join=true)
    elseif isfile(datas)
        d = [datas]
    end
    for (k, di) in enumerate(d)
        @printf("Downloading file %d out of %d... ", k, length(d))
        try
            fl = download_file("file://" * di, destionation_dir; force_download)
            preproc_data(fl, destionation_dir)
        catch e
            println(" ... error: $e")
        end
        println("done.")
    end
    nothing
end
