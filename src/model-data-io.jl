import Libz

const DEFAULT_DOWNLOAD_ROOT_URL::String = "https://www.star.bris.ac.uk/fergus/xspec/data/"
const SPECTRAL_FITTING_STORAGE_PATH = joinpath(homedir(), ".julia", "spectral_fitting_data")
const ALL_STORAGE_PATHS = String[SPECTRAL_FITTING_STORAGE_PATH]

@enum CompressionFormat::Int _CompressedGzip _NoCompression

function _extension(fmt::CompressionFormat)
    if fmt == _CompressedGzip
        return ".gz"
    elseif _NoCompression
        return ""
    else
        error("Unknown compression format $(fmt)")
    end
end

struct ModelDataInfo
    remote_path::String
    local_path::String
    compression::CompressionFormat

    function ModelDataInfo(r::String, l::String, compression::CompressionFormat)
        if compression != _NoCompression
            @assert _check_compression(r) != _NoCompression
            @assert _check_compression(l) == _NoCompression
        end
        return new(r, l, compression)
    end
end

ModelDataInfo(remote::AbstractString, local_path::AbstractString) =
    ModelDataInfo(remote, local_path, _check_compression(remote))

const _model_to_data_map = Dict{Symbol, Vector{ModelDataInfo}}()
const _model_available_memoize_cache = Dict{Symbol, Bool}()

function _check_compression(path::AbstractString)
    for f in (_CompressedGzip,)
        if endswith(path, _extension(f))
            return f
        end
    end
    return _NoCompression
end

function _trim_compression_filename(path::AbstractString)
    format = _check_compression(path)
    if format != _NoCompression
        ext = _extension(format)
        return path[1:(end - length(ext))]
    end
    return path
end

"""
    SpectralFitting.download_all_model_data()

Downloads all model data for the models currently registered with [`SpectralFitting.register_model_data`](@ref).
Calls [`SpectralFitting.download_model_data`](@ref) to perform the download.
"""
function download_all_model_data(; verbose = true)
    for s in keys(_model_to_data_map)
        download_model_data(s; verbose = verbose)
    end
    return
end

"""
    SpectralFitting.register_model_data(M::Type{<:AbstractSpectralModel}, model_data::ModelDataInfo...)
    SpectralFitting.register_model_data(M::Type{<:AbstractSpectralModel}, remote_and_local::Tuple{String,String}...)
    SpectralFitting.register_model_data(M::Type{<:AbstractSpectralModel}, filenames::String...)
    SpectralFitting.register_model_data(s::Symbol, filenames::String...)

Register `filenames` as model data associated with the model given by type `M`
or symbol `s`.  This function does not download any files, but rather adds the
relevant filenames to a lookup which
[`SpectralFitting.download_model_data`](@ref) consults when invoked, and
consequently model data is only downloaded when needed.

!!! note
    It is good practice to use this method immediately after defining a new
    model with [`@xspecmodel`](@ref) to register any required datafiles from
    the HEASoft source code, and therefore keep relevant information together.

# Example

```julia
# by type
register_model_data(XS_Laor, "ari.mod")
# by symbol
register_model_data(:XS_KyrLine, "KBHline01.fits")
```
"""
function register_model_data(s, filenames::String...; root = SPECTRAL_FITTING_STORAGE_PATH)
    model_data = map(filenames) do fname
        ModelDataInfo(fname, joinpath(root, _trim_compression_filename(fname)))
    end
    return register_model_data(_translate_model_name(s), model_data...)
end
function register_model_data(
        s,
        remote_and_local::Tuple{String, String}...;
        root = SPECTRAL_FITTING_STORAGE_PATH,
    )
    model_data = map(remote_and_local) do entry
        ModelDataInfo(entry[1], joinpath(root, entry[2]))
    end
    return register_model_data(_translate_model_name(s), model_data...)
end
function register_model_data(s::Symbol, model_data::ModelDataInfo...)
    if s in keys(_model_to_data_map)
        push!(_model_to_data_map[s], model_data...)
    else
        _model_to_data_map[s] = collect(model_data)
    end
    return
end

_translate_model_name(s::Symbol) = s
_translate_model_name(M::Type) = Base.typename(M).name

_is_model_data_downloaded(M::Type{<:AbstractSpectralModel}) =
    _is_model_data_downloaded(Base.typename(M).name)
function _is_model_data_downloaded(s::Symbol)::Bool
    if get(_model_available_memoize_cache, s, false)
        return true
    else
        return _check_model_files_and_add_to_cache(s)
    end
end

function _check_model_files_and_add_to_cache(s::Symbol)
    # get the filenames we need
    data = get(_model_to_data_map, s, nothing)
    if !isnothing(data)
        if !all(i -> ispath(i.local_path), data)
            return false
        end
    end
    # cache true
    _model_available_memoize_cache[s] = true
    return true
end

function _check_model_directory_present()
    if !ispath(SPECTRAL_FITTING_STORAGE_PATH)
        @warn "No model data directory found. Use `SpectralFitting.download_all_model_data()` to populate."
    end
    return
end

function _download_from_archive(
        src, dest;
        progress = true,
        io::IO = Core.stdout,
        model_source_url = DEFAULT_DOWNLOAD_ROOT_URL,
    )
    url = if contains(src, "://")
        src
    else
        "$model_source_url/$src"
    end
    pg = if progress
        bar = MiniProgressBar(header = "Downloading: $src", color = Base.info_color())
        start_progress(io, bar)
        (total, now) -> begin
            bar.max = total
            bar.current = now
            show_progress(io, bar)
        end
    else
        (_, _) -> nothing
    end
    Downloads.download(url, dest; progress = pg)
    end_progress(io, bar)
    return
end

"""
    SpectralFitting.download_model_data(model::AbstractSpectralModel; kwargs...)
    SpectralFitting.download_model_data(M::Type{<:AbstractSpectralModel}; kwargs...)
    SpectralFitting.download_model_data(s::Symbol; kwargs...)

Downloads the model data for a model specified either by `model`, type `M`, or symbol `s`. Datafiles
associated with a specific model may be registered using [`SpectralFitting.register_model_data`](@ref).
The download is currently unconfigurable, but permits slight control via a number of keyword arguments:

- `progress::Bool = true`

  Display a progress bar for the download.

- `model_source_url::String = "http://www.star.bris.ac.uk/fbaker/XSPEC-model-data"`

  The source URL used to download the model data.

All standard XSPEC spectral model data is currently being hosted on the University of Bristol
astrophysics servers, and should be persistently available to anyone.
"""
download_model_data(M::Type{<:AbstractSpectralModel}; kwargs...) =
    download_model_data(Base.typename(M).name; kwargs...)
download_model_data(::M; kwargs...) where {M <: AbstractSpectralModel} =
    download_model_data(Base.typename(M).name; kwargs...)
function download_model_data(
        s::Symbol;
        root = SPECTRAL_FITTING_STORAGE_PATH,
        verbose = true,
        kwargs...,
    )
    _infolog(s) = verbose && @info(s)
    if _is_model_data_downloaded(s)
        _infolog("Model data for $(s) is already downloaded.")
        return nothing
    end


    _infolog("Checking model data for $s:")
    for src in _model_to_data_map[s]
        dest = src.local_path
        mkdir_if_not_exists(dirname(dest))

        if !ispath(dest)
            mkdir_if_not_exists(dirname(dest))
            _download_from_archive(src.remote_path, dest; kwargs...)
            _infolog("$(src.local_path) downloaded")
            if src.compression != _NoCompression
                _inflate_data(src.compression, dest)
                _infolog("$(src.local_path) decompressed")
            end
        end
    end
    _infolog("All requisite model data for $s downloaded.")
    return nothing
end

function _inflate_data(compression::CompressionFormat, dest)
    if compression != _CompressedGzip
        error("Unsupported compression for deflation: $(compression)")
    end
    content = read(dest) |> Libz.ZlibInflateInputStream
    write(dest, content)
    return
end

mkdir_if_not_exists(path) = !ispath(path) && mkdir(path)

function ensure_model_data(M::Type)
    if !_is_model_data_downloaded(M)
        @warn "Model data for $(Base.typename(M).name) is not present!\nRequisite model data may be fetched with `SpectralFitting.download_model_data($(Base.typename(M).name))`."
        error("Missing data.")
    end
    return
end

function load_and_unpack_model_data(M)
    data = _model_to_data_map[Base.typename(M).name]
    # load all of the data files
    contents = map(data) do f
        path = f.local_path
        load(path)
    end
    return contents
end

"""
    get_model_data_paths(M::Type)

Get the absolute file paths to the local data files.

See also: [`get_model_data`](@ref)
"""
function get_model_data_paths(M::Type)
    ensure_model_data(M)
    data = _model_to_data_map[Base.typename(M).name]
    return [f.local_path for f in data]
end

"""
    get_model_data_paths(M::Type)

Load the model data using the `Base.load` function and return the corresponding
structures.

See also: [`get_model_data_paths`](@ref)
"""
function get_model_data(M::Type)
    ensure_model_data(M)
    return load_and_unpack_model_data(M)
end
