"""
    extract_normalisation_index(model::AbstractSpectralModel)

Return the index of the normalisation parameter (`:K`) in the
model's [`parameter_vector`](@ref). Returns `nothing` for non-additive models.

Uses the existing [`SpectralFitting.normalisationfield`](@ref) and
[`SpectralFitting.parameter_names`](@ref) from the codebase.

## Example
```julia
m = PowerLaw()
idx = extract_normalisation_index(m)  # → 1
SpectralFitting.parameter_names(m)[idx]  # → :K
```
"""
function extract_normalisation_index(
    model::AbstractSpectralModel{T, Additive},
) where {T}
    norm_sym = SpectralFitting.normalisationfield(typeof(model))
    names    = SpectralFitting.parameter_names(model)
    idx      = findfirst(==(norm_sym), names)
    isnothing(idx) && error(
        "Model $(typeof(model)) is Additive but has no field $norm_sym. " *
        "This is a bug — please report it."
    )
    return idx
end

# Non-additive models (Multiplicative, Convolutional) have no normalisation
extract_normalisation_index(::AbstractSpectralModel) = nothing

"""
    has_normalisation(model::AbstractSpectralModel)

Return `true` if `model` is an [`Additive`](@ref) model with a
normalisation parameter.
"""
has_normalisation(::AbstractSpectralModel{T, Additive}) where {T} = true
has_normalisation(::AbstractSpectralModel) = false

"""
    get_normalisation_param(model::AbstractSpectralModel{<:FitParam, Additive})

Return the [`FitParam`](@ref) for the normalisation (`:K`) of an
[`Additive`](@ref) model.
"""
get_normalisation_param(model::AbstractSpectralModel{<:FitParam, Additive}) =
    SpectralFitting.normalisation(model)

export extract_normalisation_index, has_normalisation, get_normalisation_param
