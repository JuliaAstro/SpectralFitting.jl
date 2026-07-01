"""
    struct BinnedData{Tag,V,D} <: AbstractDataset

    function BinnedData(
        domain,
        codomain;
        tag = nothing,
        domain_variance = nothing,
        codomain_variance = nothing,
        user_data = nothing,
        name = "[no-name]",
    )
        @assert length(domain) == length(codomain) + 1 "Data does not look like it is binned!"
        BinnedData{tag}(domain, codomain, domain_variance, codomain_variance, name, user_data)
    end

A datastructure for contiguously binned data (see [`ContiguouslyBinned`](@ref)).
This is the binned analogue to [`InjectiveData`](@ref).

The parameters have the following meaning:
- `Tag`: an optional user-defined tag that can be used to overwrite method dispatches (default: nothing).
- `V`: the type of the domain, co-domain, and variance vectors (default: inferred from arguments).
- `D`: the optional additional user-defined data type (default: nothing).

$(FIELDS)
"""
struct BinnedData{Tag, V, D} <: AbstractDataset
    "The domain of the data as pairs of low and high bins, such that  `E_i` and
    `E_i+1` are the low and high bins of element `i` respectively."
    domain::V
    "The co-domain (objective) of the dataset."
    codomain::V
    "The variance for each bin edge. Todo: maybe remove?"
    domain_variance::Union{Nothing, V}
    "The variance for the value in each bin."
    codomain_variance::Union{Nothing, V}
    "Dataset name."
    name::String
    "Optional user specific data that can be kept along with the dataset."
    user_data::D
end

function BinnedData(
        domain,
        codomain;
        tag = nothing,
        domain_variance = nothing,
        codomain_variance = nothing,
        user_data = nothing,
        name = "[no-name]",
    )
    @assert length(domain) == length(codomain) + 1 "Data does not look like it is binned!"
    return BinnedData{tag, promote_type(typeof(domain), typeof(codomain)), typeof(user_data)}(
        domain,
        codomain,
        domain_variance,
        codomain_variance,
        name,
        user_data,
    )
end

supports(::Type{<:BinnedData}) = (ContiguouslyBinned(), OneToOne())

bin_widths(dataset::BinnedData) = diff(dataset.domain)
objective_units(::BinnedData) = u"counts / (s * keV)"
spectrum_energy(dataset::BinnedData) =
    @views (dataset.domain[1:(end - 1)] .+ dataset.domain[2:end]) ./ 2

make_model_domain(::ContiguouslyBinned, dataset::BinnedData) = dataset.domain
make_objective(::ContiguouslyBinned, dataset::BinnedData) = dataset.codomain

make_model_domain(::OneToOne, dataset::BinnedData) = dataset.domain[1:(end - 1)]
make_objective(::OneToOne, dataset::BinnedData) = dataset.codomain

make_output_domain(layout::AbstractDataLayout, dataset::BinnedData) =
    make_model_domain(layout, dataset)

function make_objective_variance(
        ::AbstractDataLayout, dataset::BinnedData{Tag, V},
    )::V where {Tag, V}
    if !isnothing(dataset.codomain_variance)
        return dataset.codomain_variance
    else
        # TODO: i dunno just something
        return ones(eltype(V), length(dataset.codomain))
    end
end

function objective_transformer(::AbstractDataLayout, dataset::BinnedData)
    _transformer!!(domain, objective) = @views objective
    _transformer!!(output, domain, objective) = (@. output = objective)
    return _transformer!!
end

make_label(dataset::BinnedData) = dataset.name

function _printinfo(io::IO, data::BinnedData)
    printstyled(io, "BinnedData", color = :cyan)
    print(io, " with ")
    printstyled(io, length(data.codomain), color = :cyan)
    println(io, " data points:")

    dmin, dmax = prettyfloat.(extrema(data.domain))
    comin, comax = prettyfloat.(extrema(data.codomain))
    descr = """  Name                  : $(data.name)
      . Domain (min/max)    : ($dmin, $dmax)
      . Codomain (min/max)  : ($comin, $comax)
    """
    print(io, descr)
    return
end

export BinnedData
