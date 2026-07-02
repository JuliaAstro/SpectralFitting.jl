"""
    abstract type AbstractDataset

Abstract type for use in fitting routines. High level representation of some underlying
data structures.

Fitting data is considered to have an *objective* and a *domain*. As
the domain may be, for example, energy bins (high and low), or
fourier frequencies (single value), the purpose of this abstraction
is to provide some facility for translating between these representations
for the models to fit with. This is done by checking that the [`AbstractDataLayout`](@ref)
of the model and data are compatible, or at least have compatible translations.

Must implement a minimal set of accessor methods. These are paired with
`objective` and `domain` parlance. Note that these functions are prefixed with
`make_*` and not `get_*` to represent that there may be allocations or work going
into the translation. Usage of these functions should be sparse in the interest of
performance.

The arrays returned by the `make_*` functions must correspond to the [`AbstractDataLayout`](@ref)
specified by the caller.

- [`make_objective_variance`](@ref)
- [`make_objective`](@ref)
- [`make_domain_variance`](@ref)
- [`make_model_domain`](@ref)
- [`make_output_domain`](@ref)

Additionally there is an objective transformer that transforms the output of the
model onto the `output` domain:

- [`objective_transformer`](@ref)

Finally, to make all of the fitting for different statistical regimes work
efficiently, datasets should inform which units are preferred to fit. They may
also give the error statistics they prefer, and a label name primarily used to
disambiguate:

- [`preferred_units`](@ref)
- [`error_statistic`](@ref)
- [`make_label`](@ref)
"""
abstract type AbstractDataset end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(data::AbstractDataset))
    buff = IOBuffer()
    _printinfo(IOContext(buff, io), data)
    s = String(take!(buff))
    print(io, encapsulate(s))
    return
end

function Base.show(io::IO, @nospecialize(data::AbstractDataset))
    print(io, "$(Base.typename(typeof(data)).name)[$(make_label(data))]")
    return
end


"""
    make_objective(layout::AbstractDataLayout, dataset::AbstractDataset)

Returns the array used as the target for model fitting. The array must
correspond to the data [`AbstractDataLayout`](@ref) specified by the `layout`
parameter.

In as far as it can be guaranteed, the memory in the returned array will not be
mutated by any fitting procedures.

Domain for this objective should be returned by [`make_model_domain`](@ref).
"""
make_objective(layout::AbstractDataLayout, dataset::AbstractDataset) =
    error("Layout $(layout) is not implemented for $(typeof(dataset))")

"""
    make_objective_variance(layout::AbstractDataLayout, dataset::AbstractDataset)

Make the variance vector associated with each objective point.
"""
make_objective_variance(layout::AbstractDataLayout, dataset::AbstractDataset) =
    error("Layout $(layout) is not implemented for $(typeof(dataset))")

"""
    make_model_domain(layout::AbstractDataLayout, dataset::AbstractDataset)

Returns the array used as the domain for the modelling. This is paired with
[`make_domain_variance`](@ref)
"""
make_model_domain(layout::AbstractDataLayout, dataset::AbstractDataset) =
    error("Layout $(layout) is not implemented for $(typeof(dataset))")

"""
    make_domain_variance(layout::AbstractDataLayout, dataset::AbstractDataset)

Make the variance vector associated with the domain.
"""
make_domain_variance(layout::AbstractDataLayout, dataset::AbstractDataset) =
    error("Layout $(layout) is not implemented for $(typeof(dataset))")

"""
    make_output_domain(layout::AbstractDataLayout, dataset::AbstractDataset)

Returns the array used as the output domain. That is, in cases where the model
input and output map to different domains, the input domain is said to be the
model domain, the input domain is said to be the model domain.

The distinction is mainly used for the purposes of simulating data and for
visualising data.
"""
make_output_domain(layout::AbstractDataLayout, dataset::AbstractDataset) =
    error("Layout $(layout) is not implemented for $(typeof(dataset))")

"""
    objective_transformer(layout::AbstractDataLayout, dataset::AbstractDataset)

Used to transform the output of the model onto the output domain. For spectral
fitting, this is the method used to do response folding and bin masking.

If none provided, uses the `_DEFAULT_TRANSFORMER`
```julia
function _DEFAULT_TRANSFORMER()
    function _transformer!!(domain, objective)
        objective
    end
    function _transformer!!(output, domain, objective)
        @. output = objective
    end
    _transformer!!
end
```
"""
function objective_transformer(layout::AbstractDataLayout, dataset::AbstractDataset)
    @warn "Using default objective transformer!"
    return _DEFAULT_TRANSFORMER()
end

function _DEFAULT_TRANSFORMER()
    _transformer!!(domain, objective) = objective
    _transformer!!(output, domain, objective) = (@. output = objective)
    return _transformer!!
end

"""
    preferred_units(::Type{<:AbstractDataset}, s::AbstractStatistic)
    preferred_units(x, s::AbstractStatistic)

Get the preferred units that a given dataset would use to fit the
[`AbstractStatistic`](@ref) in. For example, for [`ChiSquared`](@ref), the units
of the model may be a rate, however for [`Cash`](@ref) the preferred units might
be counts.

Returning `nothing` from this function implies there is no unit preference.

If undefined for a derived type, returns `nothing`.

See also [`support_units`](@ref).
"""
preferred_units(::T, s::AbstractStatistic) where {T <: AbstractDataset} =
    preferred_units(T, s)
preferred_units(::Type{<:AbstractDataset}, s::AbstractStatistic) = nothing

"""
    error_statistic(::AbstractDataset)

Should return an [`ErrorStatistics`](@ref) describing which error statistic this
data uses.

If undefined for a derived type, returns `ErrorStatistics.Unknown`.
"""
error_statistic(::AbstractDataset) = ErrorStatistics.Unknown

"""
    make_label(d::AbstractDataset)

Return a string that gives a descriptive label for this dataset.
"""
make_label(d::AbstractDataset) = "$(Base.typename(typeof(d)).name)"

"""
Must support the same API, but may also have some query methods for specific internals.
"""
abstract type AbstractMultiDataset <: AbstractDataset end

get_datasets(data::AbstractMultiDataset) = data.datasets

function _combine_all(f, data::AbstractMultiDataset, args)
    return reduce(vcat, (f(args..., d) for d in get_datasets(data)))
end

function _printinfo(io::IO, data::AbstractMultiDataset)
    datasets = get_datasets(data)
    println(io, "MultiDataset[n=$(length(datasets))]")
    for d in datasets
        println(io, d)
    end
    return
end

make_objective(layout::AbstractDataLayout, data::AbstractMultiDataset) =
    _combine_all(make_objective, data, (layout,))
make_objective_variance(layout::AbstractDataLayout, data::AbstractMultiDataset) =
    _combine_all(make_objective_variance, data, (layout,))
make_model_domain(layout::AbstractDataLayout, data::AbstractMultiDataset) =
    _combine_all(make_model_domain, data, (layout,))
make_domain_variance(layout::AbstractDataLayout, data::AbstractMultiDataset) =
    _combine_all(make_domain_variance, data, (layout,))
make_output_domain(layout::AbstractDataLayout, data::AbstractMultiDataset) =
    _combine_all(make_output_domain, data, (layout,))

function objective_transformer(layout::AbstractDataLayout, data::AbstractMultiDataset)
    error("TODO")
end
function preferred_units(d::AbstractMultiDataset)
    error("TODO")
end
function error_statistic(d::AbstractMultiDataset)
    error("TODO")
end
function make_label(d::AbstractMultiDataset)
    error("TODO")
end

export AbstractDataset,
    error_statistic,
    make_domain_variance,
    make_label,
    make_model_domain,
    make_objective,
    make_objective_variance,
    make_output_domain,
    normalize!,
    objective_transformer,
    preferred_units

include("spectrum.jl")
include("response.jl")
include("spectraldata.jl")
include("ogipdataset.jl")
include("mission-specifics.jl")
