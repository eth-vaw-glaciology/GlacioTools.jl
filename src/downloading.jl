"""
    download_data(url, destination_dir;
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
function download_data(url::String, destination_dir::String;
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
    get_all_data(datas::Dict, destionation_dir::String;
                    force_download=false)
Downloads all files in a dictionary using `download_data`.

Input:
- datas -- dictionary of files that need to be downloaded
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
                fl = download_data(d, destionation_dir; force_download)
                preproc_data(fl, destionation_dir)
            catch e
                println(" ... error: $e")
            end
        else
            for dd in d
                try
                    fl = download_data(dd, destionation_dir; force_download)
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
