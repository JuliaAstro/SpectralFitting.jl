"""
    domain_union(datasets::AbstractDataset...; tolerance=1e-4)

Return the union of model domains from multiple datasets, merging
bin edges that are within `tolerance` of each other.

Useful when datasets from different instruments have slightly offset
energy grids that should be treated as identical bins.

## Example
```julia
union_dom = domain_union(data1, data2, data3; tolerance=1e-3)
```
"""
function domain_union(
    datasets::AbstractDataset...;
    tolerance::Float64 = 1e-4,
)
    length(datasets) >= 2 ||
        throw(ArgumentError("domain_union requires at least 2 datasets"))

    layout = ContiguouslyBinned()
    all_edges = Float64[]

    for d in datasets
        # Validate dataset supports ContiguouslyBinned
        if !(ContiguouslyBinned() in SpectralFitting.supports(typeof(d)))
            throw(ArgumentError(
                "$(typeof(d)) does not support ContiguouslyBinned layout"
            ))
        end
        append!(all_edges, make_model_domain(layout, d))
    end

    sort!(all_edges)

    # Merge edges within tolerance
    merged = Float64[]
    for e in all_edges
        if isempty(merged) || abs(e - last(merged)) > tolerance
            push!(merged, e)
        end
    end

    return merged
end

"""
    domains_are_compatible(datasets::AbstractDataset...; tolerance=1e-4)

Return `true` if all datasets share the same energy domain up to `tolerance`.

Used to verify that the single-evaluation optimisation in
[`multi_dataset_problem`](@ref) is valid.
"""
function domains_are_compatible(
    datasets::AbstractDataset...;
    tolerance::Float64 = 1e-4,
)
    layout = ContiguouslyBinned()
    domains = [make_model_domain(layout, d) for d in datasets]
    ref = first(domains)
    all(d -> length(d) == length(ref) &&
             all(abs.(d .- ref) .< tolerance), domains)
end

export domain_union, domains_are_compatible
