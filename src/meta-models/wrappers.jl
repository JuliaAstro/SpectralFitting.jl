"""
    abstract type AbstractModelWrapper{M,T,K} <: AbstractSpectralModel{T,K} end

Used to implement wrapper models that take existing models as their argument and
modify their behaviour.
"""
abstract type AbstractModelWrapper{M <: AbstractSpectralModel, T, K} <: AbstractSpectralModel{T, K} end

function remake_with_parameters(m::AbstractModelWrapper, parameters::Tuple)
    return Base.typename(typeof(m)).wrapper(remake_with_parameters(backing_model(m), parameters))
end

is_composite(::Type{<:AbstractModelWrapper{M}}) where {M} = is_composite(M)

function _inner_invokemodel!(output, domain, model::AbstractModelWrapper)
    return _inner_invokemodel!(output, domain, backing_model(model))
end

# tie in dispatches
_all_parameters_with_symbols(model::AbstractModelWrapper) =
    _all_parameters_with_symbols(backing_model(model))

backing_model(model::AbstractModelWrapper) = getfield(model, :model)

closure_and_parameter(model::AbstractModelWrapper) =
    closure_and_parameter(backing_model(model))

implementation(::Type{<:AbstractModelWrapper{M}}) where {M} = implementation(M)

parameter_count(model::AbstractModelWrapper) = parameter_count(backing_model(model))

parameter_names(::Type{<:AbstractModelWrapper{M}}) where {M} = parameter_names(M)

parameter_vector(model::AbstractModelWrapper) = parameter_vector(backing_model(model))

objective_cache_count(m::AbstractModelWrapper) = objective_cache_count(backing_model(m))

Base.copy(m::AbstractModelWrapper) = Base.typename(typeof(m)).wrapper(backing_model(m))

_wrapper_name(m::AbstractModelWrapper) = "$(Base.typename(typeof(m)).name)"
_model_name(
    m::AbstractModelWrapper,
) = "$(_wrapper_name(m))[$(_model_name(backing_model(m)))]"

_printinfo(io::IO, m::AbstractModelWrapper; kwargs...) =
    _printinfo(io, backing_model(m); kwargs...)
_printinfo(io::IO, m::AbstractModelWrapper{<:CompositeModel}; kwargs...) =
    _printinfo(io, backing_model(m); name = _wrapper_name(m), kwargs...)

remake_with_model(m::AbstractModelWrapper, model::AbstractSpectralModel) =
    Base.typename(typeof(m)).wrapper(model)

normalisation(model::AbstractModelWrapper{M, T, Additive}) where {M, T} =
    normalisation(backing_model(model))

function Base.propertynames(m::AbstractModelWrapper)
    return Base.propertynames(backing_model(m))
end

function Base.getproperty(m::AbstractModelWrapper{<:Any, <:FitParam}, symb::Symbol)
    params = Base.propertynames(backing_model(m))
    if symb in params
        return getproperty(backing_model(m), symb)
    else
        return getfield(m, symb)
    end
end

export AbstractModelWrapper
