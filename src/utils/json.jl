using CodecBzip2
using CodecZlib
using JSON

"""
    load_json(filename::AbstractString)
Load JSON data from file `filename`.
"""
function load_json(filename::AbstractString)::Dict{String,Any}
    if endswith(filename, ".json")
        data = open(filename, "r") do io
            s = read(io, String)
            JSON.parse(s)
        end
        return data
    elseif endswith(filename, ".json.bz2")
        data = open(Bzip2DecompressorStream, filename, "r") do io
            s = read(io, String)
            JSON.parse(s)
        end
        return data
    elseif endswith(filename, ".json.gz")
        data = open(GzipDecompressorStream, filename, "r") do io
            s = read(io, String)
            JSON.parse(s)
        end
        return data
    else
        error("""
        Unsupported JSON extension: \"$(basename(filename))\".
        Supported extensions are: \".json\", \".json.bz2\" and \".json.gz\".
        """)
    end
end

"""
    save_json(filename::AbstractString, data; indent)

Save `data` into JSON file `filename`. The following formats are supported:
* uncompressed JSON `.json`
* Gzip-compressed JSON `.json.gz`
* Bzip2-compressed JSON `.json.bz2`

If the file extension does not match one of the above, an error is thrown.
"""
function save_json(filename::AbstractString, data; indent=nothing)
    if endswith(filename, ".json")
        open(filename, "w") do datafile
            JSON.print(datafile, data, indent)
        end
    elseif endswith(filename, ".json.bz2")
        open(Bzip2CompressorStream, filename, "w") do datafile
            JSON.print(datafile, data, indent)
        end
    elseif endswith(filename, ".json.gz")
        open(GzipCompressorStream, filename, "w") do datafile
            JSON.print(datafile, data, indent)
        end
    else
        error("""
        Unsupported format for data file $(basename(filename)).
        Supported formats are: \".json\", \".json.bz2\" and \".json.gz\".
        """)
    end
    return nothing
end
