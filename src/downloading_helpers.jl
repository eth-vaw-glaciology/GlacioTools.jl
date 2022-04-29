"""
    download_file(url, destination_dir;
                    force_download=false,
                    filename=nothing)

Download a file, if it has not been downloaded already.
For password protected access use the `~/.netrc` file to store passwords, see
https://everything.curl.dev/usingcurl/netrc .

For downloading files on the local file system prefix them with `file://`
as you would to see them in a browser.

# Input
- url -- url for download
- destination_dir -- path of the directory to store the download

# Optional keyword args
- force_download -- force the download, even if file is present
- filename -- rename file

# Output
- destination_file -- file name (including relative path)

Example

    GlacioTools.download_file("https://raw.githubusercontent.com/pohlan/GlacioTools.jl/main/Project.toml?token=GHSAT0AAAAAABLZTWKNEUY4ZZFQ2QKU4OUQYSNSIRA", "/tmp")
"""
function download_file(url::String, destination_dir::String;
                        force_download=false,
                        filename=nothing)
    filename_in_url = split(basename(url),'?')[1] # the '?' separates query parameters, strip that too.
                                                  # TODO there might be more special chars
    mkpath(destination_dir)
    destination_file = if filename === nothing
        joinpath(destination_dir, filename_in_url)
    else
        joinpath(destination_dir, filename)
    end
    if (isdir(destination_file) || isdir(splitext(destination_file)[1]) || isfile(destination_file)) && !force_download
        # do nothing
        # print(" ... already downloaded ... ")
        return destination_file
    elseif isdir(destination_file) || isfile(destination_file)
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

TODO: .gz and .tar files?
"""
function preproc_data(fl, destination_dir)
    if splitext(fl)[2]==".zip" && !isdir(splitext(fl)[1])
        run(`unzip -ou $fl -d $destination_dir`)
        run(`rm $fl`)
    elseif splitext(fl)[2]==".gz"
        @assert splitext(splitext(fl)[1])[2]==".tar"
        run(`tar xzf $fl`)
    elseif splitext(fl)[2]==".tar"
        run(`tar xf $fl`)
    end
end

"""
    get_all_data(datas, destination_dir::String;
                    force_download=false)
Downloads all files in a directory or dictionary.

# Input
- datas -- collection of files that need to be downloaded; can be a dictionary, a folder or just the path of a single file.
- destination_dir -- path of the directory to store the download

# Optional keyword args
- force_download -- force the download, even if file is present

Example: see `fetch_Antarctica` in file `Antarctica.jl`
"""
function get_all_data(datas::Dict, destination_dir::String;
                        force_download=false)
    for (k,d) in datas
        print("Downloading $k... ")
        if d isa AbstractString
            try
                fl = download_file(d, destination_dir; force_download)
                preproc_data(fl, destination_dir)
            catch e
                println("\n error: $e")
            end
        else
            for dd in d
                try
                    fl = download_file(dd, destination_dir; force_download)
                    preproc_data(fl, destination_dir)
                catch e
                    println("\n error: $e")
                end
            end
        end
        println("done.")
    end
    nothing
end
function get_all_data(datas::String, destination_dir::String;
                        force_download=false)
    if isdir(datas)
        d = "file://" .* readdir(datas, join=true)
    elseif isfile(datas)
        d = ["file://" * datas]
    else
        d = [datas]
    end
    for (k, di) in enumerate(d)
        @printf("Downloading file %d out of %d... ", k, length(d))
        try
            fl = download_file(di, destination_dir; force_download)
            preproc_data(fl, destination_dir)
        catch e
            println("\n error: $e")
        end
        println("done.")
    end
    nothing
end
