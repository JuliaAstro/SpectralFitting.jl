# TODO: could be Response or Redistribution : how do we track this?
mutable struct ResponseMatrix{T}
    matrix::SparseMatrixCSC{T, Int}
    channels::Vector{Int}
    channel_bins_low::Vector{T}
    channel_bins_high::Vector{T}
    bins_low::Vector{T}
    bins_high::Vector{T}
end

"""
    response_energy(response::ResponseMatrix)

Get the contiguously binned energy corresponding to the *input domain* of the
response matrix. This is equivalent to the model domain.
"""
function response_energy(resp::ResponseMatrix{T}) where {T}
    E = zeros(T, length(resp.bins_low) + 1)
    E[1:(end - 1)] .= resp.bins_low
    E[end] = resp.bins_high[end]
    return E
end

"""
    folded_energy(response::ResponseMatrix)

Get the contiguously binned energy corresponding to the *output (folded) domain*
of the response matrix. That is, the channel energies as used by the spectrum.
"""
function folded_energy(resp::ResponseMatrix{T}) where {T}
    E = zeros(T, length(resp.channel_bins_low) + 1)
    E[1:(end - 1)] .= resp.channel_bins_low
    E[end] = resp.channel_bins_high[end]
    return E
end


function regroup!(resp::ResponseMatrix{T}, grouping) where {T}
    itt = GroupingIterator(grouping)
    new_matrix = zeros(T, (length(itt), size(resp.matrix, 2)))
    for grp in itt
        slice = @views resp.matrix[grp[2]:grp[3], :]
        new_matrix[grp[1], :] .= sum(slice, dims = 1) |> vec

        resp.channels[grp[1]] = grp[1] # resp.channels[grp[2]]

        resp.channel_bins_low[grp[1]] = resp.channel_bins_low[grp[2]]
        resp.channel_bins_high[grp[1]] = resp.channel_bins_high[grp[3]]
    end
    resp.matrix = sparse(new_matrix)
    return resize!(resp, length(itt))
end

function Base.resize!(response::ResponseMatrix, n)
    if size(response.matrix, 1) != n
        response.matrix = response.matrix[1:n, :]
    end
    resize!(response.channels, n)
    resize!(response.channel_bins_low, n)
    return resize!(response.channel_bins_high, n)
end

# todo: currently unused; do we need this?
function normalise_rows!(matrix)
    # returns weights
    return @views map(1:size(matrix, 2)) do
        w = sum(matrix[i, :])
        if ΣR > 0.0
            matrix[i, :] .= matrix[i, :] / w
        end
        w
    end
end

function Base.show(
        io::IO,
        ::MIME{Symbol("text/plain")},
        @nospecialize(rm::ResponseMatrix{T})
    ) where {T}
    nchans = length(rm.channels)
    println(io, "ResponseMatrix with $nchans channels:")
    Base.print_array(io, rm.matrix)
    return
end

function _printinfo(io, resp::ResponseMatrix{T}) where {T}
    emin, emax = prettyfloat(minimum(resp.bins_low)), prettyfloat(maximum(resp.bins_high))
    c_emin = prettyfloat(minimum(resp.channel_bins_low))
    c_emax = prettyfloat(maximum(resp.channel_bins_high))
    ranks, files = size(resp.matrix)

    print(io, "Response Matrix")
    printstyled(io, "($ranks x $files)", color = :cyan)
    println(io, "channels:")

    descr = """  . Chn. E (min/max)    : ($c_emin, $c_emax)
      . Domain E (min/max)  : ($emin, $emax)
    """
    print(io, descr)
    return
end

function drop_channels!(response::ResponseMatrix, inds)
    deleteat!(response.channels, inds)
    deleteat!(response.channel_bins_high, inds)
    deleteat!(response.channel_bins_low, inds)
    keep = collect(i for i in 1:size(response.matrix, 1) if i ∉ inds)
    response.matrix = response.matrix[keep, :]
    return length(inds)
end

struct AncillaryResponse{T}
    bins_low::Vector{T}
    bins_high::Vector{T}
    effective_area::Vector{T}
end

regroup!(resp::AncillaryResponse, grouping) = nothing # no op

function _printinfo(io, resp::AncillaryResponse{T}) where {T}
    emin, emax = prettyfloat(minimum(resp.bins_low)), prettyfloat(maximum(resp.bins_high))
    descr = """Ancillary Response:
      . Channels            : $(length(resp.effective_area))
      . E (min/max)         : ($emin, $emax)
    """
    print(io, descr)
    return
end

function fold_ancillary(response::ResponseMatrix, ancillary::AncillaryResponse)
    return ancillary.effective_area' .* response.matrix
end

function Base.show(
        io::IO, ::MIME{Symbol("text/plain")},
        @nospecialize(resp::AncillaryResponse{T})
    ) where {T}
    _printinfo(io, resp)
    return
end

function unfold(resp::ResponseMatrix, ancillary::AncillaryResponse, data)
    mat = fold_ancillary(resp, ancillary)
    return _unfold(resp, mat, data)
end

unfold(resp::ResponseMatrix, data) = _unfold(resp, resp.matrix, data)

function _unfold(resp::ResponseMatrix, matrix::AbstractMatrix, data)
    energy_width = (resp.bins_high .- resp.bins_low)
    bin_width = (resp.channel_bins_high .- resp.channel_bins_low)

    R = matrix * energy_width
    R_vector = vec(sum(R, dims = 2))

    return @. bin_width * (data / R_vector)
end

export ResponseMatrix
