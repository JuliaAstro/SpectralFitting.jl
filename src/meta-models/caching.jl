mutable struct CacheEntry{T}
    cache::Vector{T}
    params::Vector{T}
    domain_limits::Tuple{T, T}
    size_of_element::Int

    function CacheEntry(params::AbstractVector{<:Number})
        T = eltype(params)
        cache = zeros(T, 1)
        return new{T}(cache, params, (zero(T), zero(T)), sizeof(T))
    end
end

"""
    AutoCache

Used to automatically create a cache of another model, to avoid re-evaluating
the model if the next parameters are close to the previous parameters. The
intended use is for fitting expensive models which.

## Example

```julia
model = PhotoelectricAbsorption() * AutoCache(PowerLaw())
```

In the above model, the [`PowerLaw`](@ref) component will be augmented with the
caching behaviour.
"""
struct AutoCache{M, T, K, C <: CacheEntry} <: AbstractModelWrapper{M, T, K}
    model::M
    parameter_symbols::Vector{Symbol}
    cache::C
    abstol::Float64
    enabled::Bool
    function AutoCache(
            model::AbstractSpectralModel{T, K},
            parameter_symbols::Vector{Symbol},
            cache::CacheEntry,
            abstol,
            enabled::Bool,
        ) where {T, K}
        return new{typeof(model), T, K, typeof(cache)}(
            model,
            parameter_symbols,
            cache,
            abstol,
            enabled,
        )
    end
end

function remake_with_parameters(m::AutoCache, parameters::Tuple)
    return AutoCache(
        remake_with_parameters(backing_model(m), parameters),
        m.parameter_symbols,
        m.cache,
        m.abstol,
        m.enabled,
    )
end

function Base.copy(m::AutoCache)
    return AutoCache(
        copy(m.model),
        copy(m.parameter_symbols),
        deepcopy(m.cache),
        m.abstol,
        m.enabled,
    )
end
_model_name(model::AutoCache) = "AutoCache[$(_model_name(model.model))]"

function AutoCache(
        model::AbstractSpectralModel{T, K};
        abstol = 1.0e-3,
        enabled = true,
    ) where {T, K}
    @assert !is_composite(model)
    params, syms = _all_parameters_with_symbols(model)
    cache = CacheEntry(get_value.(params))
    return AutoCache(model, last.(syms), cache, abstol, enabled)
end

function _reinterpret_dual(
        M::Type{<:AbstractSpectralModel},
        ::Type,
        v::AbstractArray,
        n::Int,
    )
    needs_resize = n > length(v)
    if needs_resize
        @warn "$(Base.typename(M).name): Growing dual buffer..."
        resize!(v, n)
    end
    return view(v, 1:n), needs_resize
end
function _reinterpret_dual(
        M::Type{<:AbstractSpectralModel},
        DualType::Type{<:ForwardDiff.Dual},
        v::AbstractArray{T},
        n::Int,
    ) where {T}
    n_elems = div(sizeof(DualType), sizeof(T)) * n
    needs_resize = n_elems > length(v)
    if needs_resize
        @warn "$(Base.typename(M).name): Growing dual buffer..."
        resize!(v, n_elems)
    end
    return reinterpret(DualType, view(v, 1:n_elems)), needs_resize
end

function _inner_invokemodel!(
        output,
        domain,
        model::AutoCache{<:AbstractSpectralModel{T, K}},
    ) where {T, K}
    if model.enabled == false
        invoke!(output, domain, backing_model(model))
    end

    start = K === Additive ? 2 : 1
    D = promote_type(eltype(domain), T)
    _new_limits = (first(domain), last(domain))

    p_syms = getfield(model, :parameter_symbols)
    cache = getfield(model, :cache)

    # promote types for dual numbers
    output_cache, out_resized =
        _reinterpret_dual(typeof(model), D, cache.cache, length(output))
    param_cache, _ = _reinterpret_dual(typeof(model), D, cache.params, length(p_syms))

    # get the potentially new parameters of the model

    same_domain = cache.domain_limits == _new_limits

    # if the parameter size has changed, need to rerun the model
    if !out_resized && cache.size_of_element == sizeof(D) && same_domain
        # ignore the normalisation, since that's applied later
        within_tolerance = all(start:length(p_syms)) do i
            new_value = getproperty(backing_model(model), p_syms[i])
            old_value = param_cache[i]
            isapprox(new_value, old_value; rtol = getfield(model, :abstol))
        end

        if within_tolerance
            @. output = output_cache
            return output
        end
    end

    cache.size_of_element = sizeof(D)
    invoke!(output_cache, domain, backing_model(model))
    # update the auto cache infos
    cache.domain_limits = _new_limits

    # update the parameter cache
    for i in start:length(p_syms)
        param_cache[i] = getproperty(backing_model(model), p_syms[i])
    end
    # set the output
    @. output = output_cache
    return output
end

export AutoCache
