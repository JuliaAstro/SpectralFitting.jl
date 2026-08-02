"""
    multi_dataset_problem(model, datasets...; bind_norms=false)

Construct a [`FittingProblem`](@ref) for fitting a single model to multiple
datasets simultaneously, where only the normalisation (`:K`) differs between
datasets and all other parameters are shared.

Creates N copies of `model` (one per dataset), binds all shared parameters
using [`bindall!`](@ref), and returns a ready-to-fit [`FittingProblem`](@ref).

## Parameters
- `model`: An [`AbstractSpectralModel`](@ref) — copied N times
- `datasets`: Two or more [`AbstractDataset`](@ref) instances
- `bind_norms`: If `true`, also bind normalisations across datasets (default: `false`)

## Example
```julia
model = PowerLaw()
prob  = multi_dataset_problem(model, data1, data2, data3)
result = fit(prob, LevenbergMarquadt())
```
"""
function multi_dataset_problem(
    model::M,
    datasets::Vararg{<:AbstractDataset, N};
    bind_norms::Bool = false,
) where {M <: AbstractSpectralModel, N}
    N >= 2 || throw(ArgumentError("Need at least 2 datasets, got $N"))

    # Type-stable: ntuple with Val(N) gives Tuple, not Vector
    models = ntuple(_ -> deepcopy(model), Val(N))

    # Type-stable pair construction
    prob = FittingProblem(map(=>, models, datasets)...)

    # Get shared parameter names (all except normalisation)
    all_names  = SpectralFitting.parameter_names(model)
    norm_sym   = has_normalisation(model) ?
                 SpectralFitting.normalisationfield(typeof(model)) : nothing

    shared_names = filter(s -> s !== norm_sym, all_names)

    # Bind shared params across all N models
    for sym in shared_names
        bindall!(prob, sym)
    end

    if bind_norms && !isnothing(norm_sym)
        bindall!(prob, norm_sym)
    end

    return prob
end

export multi_dataset_problem
