"""
    AsConvolution(model; domain = collect(range(0, 2, 100)))

Turn an additive model into a convolutional model.

## Example

```
convolution_model = AsConvolution(GaussianLine())
```

The above model will now convolve the [`GaussianLine`](@ref) model onto whatever
it is applied to.

The `domain` keyword can be used to pass the domain on which to evaluate the
wrapped model before convolution. The domain is such that `x = 1` corresponds to
no domain shift in the convolution (e.g., if the domain is energy, this would be
the rest energy of a line).
"""
struct AsConvolution{M, T, V, P} <: AbstractModelWrapper{M, T, Convolutional}
    model::M
    # the domain on which we evaluate this model
    domain::V
    # an additional output cache
    cache::NTuple{2, Vector{P}}
    function AsConvolution(
            model::AbstractSpectralModel{T},
            domain::V,
            cache::NTuple{2, Vector{P}},
        ) where {T, V, P}
        @assert !is_composite(model)
        return new{typeof(model), T, V, P}(model, domain, cache)
    end
end

is_composite(::Type{<:AsConvolution}) = false

"""
    Base.copy(m::AsConvolution)

Creates a copy of an [`AsConvolution`](@ref) wrapped model. Will make a
`deepcopy` of the cache to eliminate possible thread contention, but does not
copy the domain.
"""
Base.copy(m::AsConvolution) =
    AsConvolution(copy(backing_model(m)), m.domain, deepcopy(m.cache))

# ignore K
function parameter_names(
        ::Type{<:AsConvolution{M}},
    ) where {M <: AbstractSpectralModel{T, Additive}} where {T}
    return parameter_names(M)[2:end]
end
function parameter_count(m::AsConvolution{<:AbstractSpectralModel{T, Additive}}) where {T}
    return parameter_count(backing_model(m)) - 1
end
function parameter_vector(
        model::AsConvolution{<:AbstractSpectralModel{T, Additive}},
    ) where {T}
    return parameter_vector(backing_model(model))[2:end]
end

# tie in dispatches
function _all_parameters_with_symbols(
        model::AsConvolution{<:AbstractSpectralModel{T, Additive}},
    ) where {T}
    ps, syms = _all_parameters_with_symbols(backing_model(model))
    return ps[2:end], syms[2:end]
end

function remake_with_parameters(
        model::AsConvolution{<:AbstractSpectralModel{T, K}},
        params::Tuple,
    ) where {T, K}
    _params = if K <: Additive
        # Need an additional parameter for the normalisation term
        (one(eltype(params)), params...)
    else
        params
    end
    return AsConvolution(remake_with_parameters(model.model, _params), model.domain, model.cache)
end

_model_name(model::AsConvolution) = "AsConvolution[$(_model_name(model.model))]"

function AsConvolution(
        model::AbstractSpectralModel{T};
        domain = collect(range(0, 2, 100)),
    ) where {T}
    output = collect(invokemodel(domain, model))
    return AsConvolution(model, domain, (output, deepcopy(output)))
end

function _inner_invokemodel!(output, domain, model::AsConvolution{M, T}) where {M, T}
    D = promote_type(eltype(domain), T)
    model_output, _ =
        _reinterpret_dual(typeof(model), D, model.cache[1], length(model.domain) - 1)
    convolution_cache, _ =
        _reinterpret_dual(typeof(model), D, model.cache[2], length(output))

    # invoke the child model
    invoke!(model_output, model.domain, model.model)

    # do the convolution
    convolve!(convolution_cache, output, domain, model_output, model.domain)

    # overwrite the output
    @. output = convolution_cache
    return
end

export AsConvolution
