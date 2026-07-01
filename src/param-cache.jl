# these should work for both single models and composite models
# so that the parameters can all be allocated in one go
# for multi model: give a view to each cache
struct ParameterCache{M <: AbstractArray, V, T <: Number}
    free_mask::M # bit vector or a view into one
    parameters::V
    frozen_values::Vector{T}
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(pc::ParameterCache))
    println(io, "ParameterCache")
    println(io, " . Free Mask     : ", pc.free_mask)
    println(io, " . Parameters    : ", _get_parameters(pc, 0))
    println(io, " . Frozen Values : ", pc.frozen_values)
    return
end

function _make_free_mask(params::AbstractArray{<:FitParam})
    free_mask = BitVector(undef, length(params))
    for (i, p) in enumerate(params)
        free_mask[i] = isfree(p)
    end
    return free_mask
end

function ParameterCache(params::AbstractArray{<:FitParam})
    free_mask = _make_free_mask(params)
    frozen = params[.!free_mask]
    return ParameterCache(free_mask, map(get_value, params), map(get_value, frozen))
end

_get_parameters(cache::ParameterCache, params) = cache.parameters
_get_parameters(cache::ParameterCache{M, V}, params) where {M <: AbstractArray, V <: DiffCache} =
    get_tmp(cache.parameters, params)

function update_free_parameters!(cache::ParameterCache, params)
    @assert count(cache.free_mask) == length(params) "$(count(cache.free_mask)) != $(length(params))"

    _params = _get_parameters(cache, zero(eltype(params)))
    # copy over the frozen parameters as well
    j::Int = 1
    k::Int = 1
    for (i, mask) in enumerate(cache.free_mask)
        if !mask
            _params[i] = cache.frozen_values[j]
            j += 1
        else
            _params[i] = params[k]
            k += 1
        end
    end

    return _params
end

function get_free_parameters(cache::ParameterCache)
    return _get_parameters(cache, 0)[cache.free_mask]
end

export ParameterCache, update_free_parameters!
