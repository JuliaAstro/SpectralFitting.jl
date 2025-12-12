module OGIP

import SpectralFitting
import SpectralFitting: SpectralUnits

import FITSFiles

using SparseArrays

struct MissingHeader <: Exception
    header::String
end
Base.showerror(io::IO, e::MissingHeader) = print(io, "Header: '$(e.header)' is not defined")

# structs
struct RMFHeader
    first_channel::Int
    num_channels::Int
end

struct RMFMatrix{V,T,M}
    f_chan::V
    n_chan::V
    bins_low::Vector{T}
    bins_high::Vector{T}
    matrix::M
    header::RMFHeader
end

struct RMFChannels{T}
    channels::Vector{Int}
    bins_low::Vector{T}
    bins_high::Vector{T}
end

function _parse_any(::Type{T}, @nospecialize(value::V))::T where {T,V}
    if V <: AbstractString
        parse(T, value)
    else
        convert(T, value)
    end
end

function _string_boolean(@nospecialize(value::V))::Bool where {V}
    if V <: AbstractString
        if value == "F"
            false
        elseif value == "T"
            true
        else
            @warn("Unknown boolean string: $(value)")
            false
        end
    else
        value
    end
end

# functions
function parse_rmf_header(table::FITSFiles.HDU)
    findex = findfirst(==(:F_CHAN), keys(table.data))
    if isnothing(findex)
        throw(MissingHeader("F_CHAN"))
    end

    tlindex = "TLMIN$findex"
    first_channel = if haskey(table.cards, tlindex)
        _parse_any(Int, get(table.cards, tlindex))
    else
        @warn "No TLMIN key set in RMF header ($tl_min_key). Assuming channels start at 1."
        1
    end
    num_channels = if haskey(table.cards, "DETCHANS")
        _parse_any(Int, get(table.cards, "DETCHANS"))
    else
        @warn "DETCHANS is not set in RMF header. Infering channel count from table length."
        -1
    end
    RMFHeader(first_channel, num_channels)
end

function read_rmf_channels(table::FITSFiles.HDU, T::Type)
    channels = _parse_any.(Int, table.data.CHANNEL)
    energy_low = _parse_any.(T, table.data.E_MIN)
    energy_high = _parse_any.(T, table.data.E_MAX)
    RMFChannels(channels, energy_low, energy_high)
end

function _chan_to_vectors(chan::AbstractMatrix)
    map(eachcol(chan)) do column
        i = findfirst(==(0), column)
        # if no zeroes, return full column
        if isnothing(i)
            column
        else
            # exclude the zero from our selection
            i = i > 1 ? i - 1 : i
            column[1:i]
        end
    end
end

function _translate_channel_array(channel)
    if channel isa AbstractMatrix
        # transpose the matrix, since FITSFiles reads things in Julia-style
        _chan_to_vectors(transpose(channel))
    elseif eltype(channel) <: AbstractVector
        # reorder for the same reason as above transpose
        channel
    else
        map(i -> [i], channel)
    end
end

function _adapt_matrix_type(T::Type, mat::M) where {M}
    if eltype(M) <: AbstractVector
        map(row -> convert.(T, row), mat)
    elseif M <: AbstractMatrix
        map(row -> convert.(T, row), eachcol(mat))
    else
        # handle simple vector format where each energy bin maps to a single
        # value (e.g., for responses produced using ftflx2xsp)
        map(val -> [convert(T, val)], mat)
    end
end

function read_rmf_matrix(table::FITSFiles.HDU, header::RMFHeader, T::Type)
    energy_low = convert.(T, table.data.ENERG_LO)
    energy_high = convert.(T, table.data.ENERG_HI)
    f_chan_raw = table.data.F_CHAN
    n_chan_raw = table.data.N_CHAN
    matrix_raw = table.data.MATRIX

    # type stable: convert to common vector of vector format
    f_chan::Vector{Vector{Int}} = _translate_channel_array(f_chan_raw)
    n_chan::Vector{Vector{Int}} = _translate_channel_array(n_chan_raw)

    RMFMatrix(
        f_chan,
        n_chan,
        energy_low,
        energy_high,
        _adapt_matrix_type(T, matrix_raw),
        header,
    )
end

# TODO: remove me
function _read_fits_and_close(f, path)
    fits_file = FITSFiles.fits(path)
    f(fits_file)
end

function read_rmf(path::String; T::Type = Float64)
    (header, rmf, channels::RMFChannels{T}) = _read_fits_and_close(path) do fits
        rmf_i = find_extension(fits, ["RESP", "MATRIX"])
        energy_i = find_extension(fits, "EBOUND")
        hdr = parse_rmf_header(fits[rmf_i])
        _rmf = read_rmf_matrix(fits[rmf_i], hdr, T)
        _channels = read_rmf_channels(fits[energy_i], T)
        (hdr, _rmf, _channels)
    end

    _build_reponse_matrix(header, rmf, channels, T)
end

function find_extension(
    fits,
    extension::T,
) where {T<:Union{<:AbstractString,<:AbstractVector}}
    # find the correct extensions
    i::Int = 2
    for hdu in fits
        extname = get(hdu.cards, "EXTNAME", nothing)
        if isnothing(extname)
            continue
        end

        if T <: AbstractString
            if contains(extname, extension)
                return i
            end
        elseif T <: AbstractVector
            for ext in extension
                if contains(extname, ext)
                    return i
                end
            end
        end
        i += 1
    end
    return nothing
end

function read_ancillary_response(path::String; T::Type = Float64)
    (bins_low, bins_high, effective_area) = _read_fits_and_close(path) do fits
        i = find_extension(fits, "RESP")
        area::Vector{T} = convert.(T, fits[i].data.SPECRESP)
        lo::Vector{T} = convert.(T, fits[i].data.ENERG_LO)
        hi::Vector{T} = convert.(T, fits[i].data.ENERG_HI)
        (lo, hi, area)
    end
    SpectralFitting.AncillaryResponse{T}(bins_low, bins_high, effective_area)
end

function _build_reponse_matrix(
    header::RMFHeader,
    rmf::RMFMatrix,
    channels::RMFChannels,
    T::Type,
)
    R = build_response_matrix(
        rmf.f_chan,
        rmf.n_chan,
        rmf.matrix,
        header.num_channels,
        header.first_channel,
        T,
    )
    SpectralFitting.ResponseMatrix(
        R,
        channels.channels,
        channels.bins_low,
        channels.bins_high,
        rmf.bins_low,
        rmf.bins_high,
    )
end

function build_response_matrix(
    f_chan::Vector,
    n_chan::Vector,
    matrix_rows::Vector,
    num_cols::Int,
    first_channel,
    T::Type,
)
    ptrs = Int[1]
    indices = Int[]
    matrix = Float64[]

    prev = first(ptrs)

    for i in eachindex(f_chan)
        M = matrix_rows[i]

        row_len = 0
        for (f, n) in zip(f_chan[i], n_chan[i])
            if n == 0
                # advance row
                break
            end
            first = (f - first_channel) + 1
            # append all of the indices
            for j = 0:(n-1)
                push!(indices, j + first)
            end
            append!(matrix, M[(row_len+1):(row_len+n)])
            row_len += n
        end

        next = row_len + prev
        push!(ptrs, next)
        prev = next
    end

    SparseArrays.SparseMatrixCSC{T,Int}(num_cols, length(f_chan), ptrs, indices, matrix)
end

# TODO: marked for refactoring (unused)
function build_response_matrix!(
    R,
    f_chan::Vector,
    n_chan::Vector,
    matrix_rows::Vector,
    first_channel,
)
    for (i, (F, N)) in enumerate(zip(f_chan, n_chan))
        M = matrix_rows[i]
        index = 1
        for (first, len) in zip(F, N)
            if len == 0
                break
            end
            first -= first_channel
            @views R[(first+1):(first+len), i] .= M[index:(index+len-1)]
            index += len
        end
    end
end

function _get_exposure_time(header)
    if haskey(header, "EXPOSURE")
        return get(header, "EXPOSURE")
    end
    if haskey(header, "TELAPSE")
        return get(header, "TELAPSE")
    end
    # maybe time stops given
    if (haskey(header, "TSTART")) && (haskey(header, "TSTOP"))
        return get(header, "TSTOP") - get(header, "TSTART")
    end
    @warn "Cannot find or infer exposure time."
    0.0
end

function _get_stable(::Type{T}, header, name, default)::T where {T}
    get(header, name, T(default))
end

function read_spectrum(path; T::Type = Float64)
    info::SpectralFitting.Spectrum{T} = _read_fits_and_close(path) do fits
        header = fits[2].cards
        # if not set, assume not poisson errors
        is_poisson = _string_boolean(get(header, "POISSERR", false))
        # read general infos
        instrument = strip(get(header, "INSTRUME"))
        telescope = strip(get(header, "TELESCOP"))
        exposure_time = T(_get_exposure_time(header))
        background_scale = _get_stable(T, header, "BACKSCAL", one(T))
        area_scale = _get_stable(T, header, "AREASCAL", one(T))
        sys_error = _get_stable(T, header, "SYS_ERR", zero(T))

        column_names = keys(fits[2].data)

        channels::Vector{Int} = convert.(Int, fits[2].data.CHANNEL)

        quality::Vector{Int} = if :QUALITY ∈ column_names
            convert.(Int, fits[2].data.QUALITY)
        else
            zeros(Int, size(channels))
        end

        grouping::Vector{Int} = if :GROUPING ∈ column_names
            convert.(Int, fits[2].data.GROUPING)
        else
            ones(Int, size(channels))
        end

        units::SpectralUnits.RateOrCount, values::Vector{T} = if :RATE ∈ column_names
            SpectralUnits._rate(), convert.(T, fits[2].data.:RATE)
        else
            SpectralUnits._counts(), convert.(T, fits[2].data.:COUNTS)
        end

        stat, _errors = if :STAT_ERR ∈ column_names
            if is_poisson
                @warn "Both STAT_ERR column present and POISSERR flag set. Using STAT_ERR."
            end
            SpectralFitting.ErrorStatistics.Numeric, convert.(T, fits[2].data.STAT_ERR)
        elseif is_poisson
            SpectralFitting.ErrorStatistics.Poisson,
            @. T(SpectralFitting.count_error(values, 1.0))
        else
            @warn "Unknown error statistics. Setting zero for all."
            SpectralFitting.ErrorStatistics.Unknown, T[0 for _ in values]
        end

        SpectralFitting.Spectrum{T}(
            channels,
            quality,
            grouping,
            values,
            units,
            exposure_time,
            background_scale,
            area_scale,
            stat,
            _errors,
            sys_error,
            telescope,
            instrument,
        )
    end
    info
end

function read_background(path::String)
    read_spectrum(path)
end

function read_paths_from_spectrum(path::String)
    header = _read_fits_and_close(path) do fits
        # TODO: FITSFiles will here parse the whole file just so we can get at
        # the header, which seems a bit redundant...
        fits[2].cards
    end

    # extract path information from header
    possible_ext = splitext(path)[2]
    response_path = read_filename(header, "RESPFILE", path, ".rmf", ".rsp")
    ancillary_path = read_filename(header, "ANCRFILE", path, possible_ext)
    background_path = read_filename(header, "BACKFILE", path, possible_ext)
    (background_path, response_path, ancillary_path)
end

function read_filename(header, entry, parent, exts...)
    data_directory = Base.dirname(parent)
    parent_name = basename(parent)
    if haskey(header, entry)
        path::String = strip(get(header, entry))
        if path == "NONE"
            return nothing
        end
        name = find_file(data_directory, path, parent_name, exts)
        if !isnothing(name)
            return name
        end
    end
    nothing
end

function find_file(dir, name, parent, extensions)
    if length(name) == 0
        return nothing
    elseif match(r"%match%", name) !== nothing
        base = splitext(parent)[1]
        for ext in extensions
            testfile = joinpath(dir, base * ext)
            if isfile(testfile)
                return testfile
            end
        end
        @warn "Missing! Could not find file '%match%': tried $extensions"
        return nothing
    elseif match(r"^none\b", name) !== nothing
        return nothing
    end
    joinpath(dir, name)
end

end # module

using .OGIP
export OGIP

function read_fits_header(path; hdu = 2)
    OGIP._read_fits_and_close(path) do f
        # TODO: FITSFiles will here parse the whole file just so we can get at
        # the header, which seems a bit redundant...
        f[hdu].cards
    end
end

export read_fits_header
