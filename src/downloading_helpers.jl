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
                                                  # TODO there might be more special URL-chars?
    mkpath(destination_dir)
    destination_file = if filename === nothing
        joinpath(destination_dir, filename_in_url)
    else
        joinpath(destination_dir, filename)
    end
    @assert destination_dir == dirname(destination_file)
    if isfile(destination_file) && !force_download
        # do nothing
        # print(" ... already downloaded ... ")
    elseif isfile(destination_file)
        rm(destination_file)
        Downloads.download(url, destination_file)

        # if startswith(url, "file://")
        #     # just copy, without the "file://" part
        #     cp(url[8:end], destination_file, force=force_download)
        # else
        #     Downloads.download(url, destination_file)
        # end
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
                        force_download=false,
                        preproc=true)

Downloads all files in a directory or dictionary and runs `preproc_data` on them too.

# Input
- datas -- collection of files that need to be downloaded; can be a dictionary,
           a folder-path or just the path of a single file.
           --> folders are not downloaded recursively!
- destination_dir -- path of the directory to store the download

# Optional keyword args
- force_download -- force the download, even if file is present [false]
- preproc -- run `preproc_data` [true]

Example: see `fetch_Antarctica` in file `Antarctica.jl`
"""
function get_all_data(datas::Dict, destination_dir::String;
                      force_download=false,
                      preproc=true)
    for (k,d) in datas
        print("Downloading $k... ")
        if d isa AbstractString
            try
                fl = download_file(d, destination_dir; force_download)
                preproc && preproc_data(fl, destination_dir)
            catch e
                println("\n error: $e")
            end
        elseif d isa Vector
            for dd in d
                try
                    @show fl = download_file(dd, destination_dir; force_download)
                    preproc && preproc_data(fl, destination_dir)
                catch e
                    println("\n error: $e")
                end
            end
        else
            error("`datas` Dict values need to be strings or list of strings")
        end
        println("done.")
    end
    nothing
end
function get_all_data(datas::String, destination_dir::String;
                      force_download=false,
                      preproc=true)
    if isdir(datas)
        d = ["file://" * fl for fl in readdir(datas, join=true) if isfile(fl)]
    elseif isfile(datas)
        d = ["file://" * datas]
    else
        d = [datas]
    end
    for (k, di) in enumerate(d)
        @printf("Downloading file %d out of %d... ", k, length(d))
        try
            fl = download_file(di, destination_dir; force_download)
            preproc && preproc_data(fl, destination_dir)
        catch e
            println("\n error: $e")
        end
        println("done.")
    end
    nothing
end
