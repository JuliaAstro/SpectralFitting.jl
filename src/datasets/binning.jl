export construct_objective_cache, regroup!, interpolated_rebin!, interpolated_rebin

function interpolated_rebin(target_bins, input, current_bins)
    output = zeros(eltype(input), length(target_bins_high) - 1)
    return interpolated_rebin!(output, target_bins, input, current_bins)
end

function interpolated_rebin!(output, target, input, current)
    @assert length(output) == length(target) - 1
    @assert length(input) == length(current) - 1

    # ensure everything is zeroed
    @. output = 0

    # convert from counts to density
    for i in eachindex(input)
        input[i] /= current[i + 1] - current[i]
    end

    # use DataInterpolations cus no allocation :D
    interp = @views DataInterpolations.LinearInterpolation(input, current[1:(end - 1)])

    function _interp_safe(x)
        return x < current[1] || x > current[end - 1] ? zero(x) : interp(x)
    end

    for i in eachindex(output)
        bin_high, bin_low = target[i + 1], target[i]
        Δ = bin_high - bin_low
        output[i] = Δ * _interp_safe(bin_low)
    end

    # restore like a good mannered citizen
    for i in eachindex(input)
        input[i] *= current[i + 1] - current[i]
    end
    return output
end

function construct_objective_cache(
        T::Type,
        model::AbstractSpectralModel,
        domain::AbstractVector,
    )
    return construct_objective_cache(preferred_support(model), T, model, domain)
end
function construct_objective_cache(
        layout::AbstractDataLayout,
        model::AbstractSpectralModel,
        domain::AbstractVector{T},
    ) where {T}
    return construct_objective_cache(layout, T, model, domain)
end
function construct_objective_cache(
        layout::AbstractDataLayout,
        T::Type,
        model::AbstractSpectralModel,
        domain,
    )
    return construct_objective_cache(layout, T, length(domain), objective_cache_count(model))
end
function construct_objective_cache(layout::AbstractDataLayout, T::Type, N::Int, size::Int)
    n = if layout isa OneToOne
        N
    elseif layout isa ContiguouslyBinned
        N - 1
    else
        error("No domain construction for $(layout) defined.")
    end
    return zeros(T, (n, size))
end

function augmented_energy_channels(channels, other_channels, bins_low, bins_high)
    # TODO: just going to assume the channels line up
    N = length(channels)
    energies = zeros(eltype(bins_high), N + 1)
    for (i, c) in enumerate(channels)
        index = findnext(==(c), other_channels, i)
        if isnothing(index)
            error(
                "Failed to find channel in response corresponding to channel $c in spectrum.",
            )
        end
        if index > lastindex(bins_low)
            break
        end
        if i > 1 && !(energies[i] ≈ bins_low[index])
            @warn "Channel $index: misaligned $(energies[i - 1]) != $(bins_low[index])! Data may not be contiguous."
        end
        energies[i] = bins_low[index]
        energies[i + 1] = bins_high[index]
    end
    return energies
end
