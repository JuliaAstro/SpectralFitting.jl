struct DomainExtension{T}
    low::Vector{T}
    high::Vector{T}
end

function DomainExtension(T::Type)
    return DomainExtension(T[], T[])
end

"""
    FittableMultiDataset(d1, d2...)

A thin wrapper representing multiple datasets.
"""
struct FittableMultiDataset{D, T}
    d::D
    extension::Vector{DomainExtension{T}}
    function FittableMultiDataset(data::Vararg{<:AbstractDataset}; T::Type = Float64)
        return new{typeof(data), T}(data, [DomainExtension(T) for _ in data])
    end
end

function Base.getindex(multidata::FittableMultiDataset, i)
    return multidata.d[i]
end

function Base.merge(m1::FittableMultiDataset, m2::FittableMultiDataset)
    return FittableMultiDataset(m1.d..., m2.d...)
end

"""
    FittableMultiModel(m1, m2...)

A thin wrapper representing multiple models.
"""
struct FittableMultiModel{M}
    m::M
    FittableMultiModel(model::Vararg{<:AbstractSpectralModel}) = new{typeof(model)}(model)
end

function Base.getindex(multimodel::FittableMultiModel, i)
    return multimodel.m[i]
end

function _multi_constructor_wrapper(
        T::Union{<:Type{<:FittableMultiDataset}, <:Type{<:FittableMultiModel}},
        args,
    )
    if args isa T
        return args
    end
    if args isa Tuple
        if eltype(args) <: T
            return reduce(merge, args)
        else
            return T(args...)
        end
    end
    return T(args)
end

"""
    _accumulated_indices(items)

`items` is a tuple or vector of lengths `n1, n2, ...`

Returns a tuple or array with same length as items, which gives the index boundaries of
an array with size `n1 + n2 + ...`.
"""
function _accumulated_indices(items)
    total::Int = 0
    return map(items) do item
        total = item + total
        total
    end
end

"""
    ParameterTriple(model_index, component, parameter)
    ParameterTriple(info_tuple)

Used to unambiguously denote a given parameter in a fitting problem.
"""
struct ParameterTriple
    "Model index number."
    model_index::Int
    "If the model is a composite model, the symbol of the composite component."
    component::Union{Nothing, Symbol}
    "The symbol of the parameter."
    parameter::Symbol
end

function ParameterTriple(t::Tuple)
    N = length(t)
    @assert 2 <= N <= 3 "Tuple has bad length ($(length(t)))"
    if N == 2
        return ParameterTriple(t[1], nothing, t[2])
    else
        return ParameterTriple(t[1], t[2], t[3])
    end
end

function _build_triples(model::AbstractSpectralModel, model_index::Int)
    if is_composite(model)
        out = ParameterTriple[]
        for m in propertynames(model)
            for s in parameter_names(getproperty(model, m))
                push!(out, ParameterTriple(model_index, m, s))
            end
        end
    else
        out = [ParameterTriple(model_index, nothing, s) for s in parameter_names(model)]
    end
    return out
end

function _build_parameter_lookup(m::FittableMultiModel)
    lookup = Dict{ParameterTriple, Int}()
    i = 1
    for (model_index, model) in enumerate(m.m)
        if is_composite(model)
            for component in propertynames(model)
                for sym in parameter_names(getproperty(model, component))
                    lookup[ParameterTriple(model_index, component, sym)] = i
                    i += 1
                end
            end
        else
            for sym in parameter_names(model)
                lookup[ParameterTriple(model_index, nothing, sym)] = i
                i += 1
            end
        end
    end
    return lookup
end

"""
    FittingProblem

A struct representing a combination of models and dataset to be fit.
"""
struct FittingProblem{M <: FittableMultiModel, D <: FittableMultiDataset}
    model::M
    data::D
    "Given a parameter triple, returns the index of the corresponding parameter."
    lookup::Dict{ParameterTriple, Int}
    bindings::Dict{Int, Vector{Int}}

    function FittingProblem(m::FittableMultiModel, d::FittableMultiDataset)
        lookup = _build_parameter_lookup(m)
        return new{typeof(m), typeof(d)}(m, d, lookup, Dict{Int, Vector{Int}}())
    end
end

FittingProblem(pairs::Vararg{<:Pair}) = FittingProblem(first.(pairs), last.(pairs))

simplify!(prob::FittingProblem) = _simplify_bindings!(prob.bindings)

function translate_bindings(prob::FittingProblem)
    revlookup(param_index) = only([k for (k, v) in prob.lookup if v == param_index])

    IndexType = Dict{Symbol, Dict{Symbol, String}}
    model_translation = Dict{Int, IndexType}()

    for (root, targets) in prob.bindings
        main = revlookup(root)
        main_str = if isnothing(main.component)
            "Bound to Model $(main.model_index) -> $(main.parameter)"
        else
            "Bound to Model $(main.model_index) => $(main.component) -> $(main.parameter)"
        end

        for t in targets
            target = revlookup(t)

            existing = get(model_translation, target.model_index, IndexType())

            _component = if isnothing(target.component)
                :only
            else
                target.component
            end

            index = get(existing, _component, Dict{Symbol, String}())

            index[target.parameter] = main_str

            existing[_component] = index
            model_translation[target.model_index] = existing
        end
    end

    return model_translation
end

function _simplify_bindings!(bindings::Dict{Int, Vector{Int}})
    isempty(bindings) && return

    roots = sort!(collect(keys(bindings)))
    skip = Set()

    for root in roots
        root in skip && continue
        targets = bindings[root]

        for t in targets
            if haskey(bindings, t)
                push!(skip, t)
                targets = vcat(targets, bindings[t])
                delete!(bindings, t)
            end
        end

        bindings[root] = unique!(targets)
    end

    return bindings
end

function _build_parameter_mapping(prob::FittingProblem{M}) where {M <: FittableMultiModel}
    simplify!(prob)
    pvec, _ = _all_parameters_with_symbols(prob)
    total_params = length(pvec)

    reverse_bindings = Dict{Int64, Int64}()
    for i in prob.bindings
        for j in prob.bindings[i[1]]
            reverse_bindings[j] = i[1]
        end
    end

    bound_parameters = keys(reverse_bindings)
    parameter_bindings = collect(1:total_params)

    for i in bound_parameters
        parameter_bindings[i] = reverse_bindings[i]
    end

    I = (1:length(M.parameters[1].parameters)...,)
    _start = Ref{Int}(1)
    return map(I) do model_index
        start = _start[]
        model = prob.model.m[model_index]
        N = length(parameter_vector(model))
        b = parameter_bindings[start:(start + N - 1)]
        _start[] += N
        b
    end
end

function parameter_vector_symbols_and_bindings(prob::FittingProblem)
    pvec, syms = _all_parameters_with_symbols(prob)
    bindings = _build_parameter_mapping(prob)
    return pvec, syms, bindings
end

function _all_parameters_with_symbols(prob::FittingProblem)
    # TODO: I don't like that this has a different return type
    symbol_mappings = Pair{Symbol, Symbol}[]
    params = map(prob.model.m) do m
        ps, syms = _all_parameters_with_symbols(m)
        for pair in syms
            push!(symbol_mappings, pair)
        end
        ps
    end
    return reduce(vcat, params), symbol_mappings
end

function FittingProblem(m, d)
    _m = _multi_constructor_wrapper(FittableMultiModel, m)
    _d = _multi_constructor_wrapper(FittableMultiDataset, d)
    return FittingProblem(_m, _d)
end

function model_count(prob::FittingProblem)
    return length(prob.model.m)
end

function data_count(prob::FittingProblem)
    return length(prob.data.d)
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(prob::FittingProblem))
    buff = IOBuffer()
    buff_c = IOContext(buff, io)

    println(buff_c, "FittingProblem:")
    println(buff_c, "  . Models     : $(length(prob.model.m))")
    println(buff_c, "  . Datasets   : $(length(prob.data.d))")
    println(buff_c, "  Parameter Summary:")

    factor = div(length(prob.data.d), length(prob.model.m))

    total = 0
    free = 0
    frozen = 0
    for m in prob.model.m
        for p in parameter_vector(m)
            total += factor
            if isfrozen(p)
                frozen += factor
            else
                free += factor
            end
        end
    end

    bound = length(unique!(reduce(vcat, values(prob.bindings), init = Int[])))

    free = free - bound

    println(buff_c, "  . Total      : $(total)")

    print(buff_c, "  . ")
    printstyled(buff_c, "Frozen", color = :cyan)
    println(buff_c, "     : $(frozen)")

    print(buff_c, "  . ")
    printstyled(buff_c, "Bound", color = :cyan)
    println(buff_c, "      : $(bound)")

    print(buff_c, "  . ")
    printstyled(buff_c, "Free", color = :cyan)
    println(buff_c, "       : $(free)")

    print(io, encapsulate(String(take!(buff))))
    return
end


"""
    details(prob::FittingProblem)

Show details about the fitting problem, including the specific model parameters
that are bound together.
"""
function details(prob::FittingProblem; color = true)
    buff = IOBuffer()
    buff_ctx = IOContext(buff, :color => color)

    t_binds = translate_bindings(prob)

    print(buff_ctx, "Models:")
    for (i, m) in enumerate(prob.model.m)
        buf = IOBuffer()
        buf_ctx = IOContext(buf, buff_ctx)

        println(buf_ctx)
        printstyled(buf_ctx, "Model $i", color = :yellow)
        print(buf_ctx, ": ")


        _bindings = if is_composite(m)
            get(t_binds, i, nothing)
        else
            t_model = get(t_binds, i, nothing)
            if isnothing(t_model)
                nothing
            else
                t_model[:only]
            end
        end

        _printinfo(buf_ctx, m; bindings = _bindings)
        print(buff_ctx, indent(String(take!(buf)), 2))
    end

    print(encapsulate(String(take!(buff))))
    return
end

_sort_binding!(binding) = sort!(binding, by = i -> i[1])

function _construct_bound_mapping(bindings, parameter_count)
    remove = Int[]

    parameter_mapping = map((1:length(parameter_count)...,)) do i
        # indices of the full parameter array for the ith model
        collect(range(_get_range(parameter_count, i)...))
    end

    for binding in bindings
        # the model we use as a reference to bind *to*
        reference = binding[1]
        for b in @views binding[2:end]
            # get the number of the parameter that we are binding
            parameter_number = parameter_mapping[b[1]][b[2]]
            parameter_mapping[b[1]][b[2]] = reference[2]

            # mark for removal: find the parameter index in the global array
            N = if b[1] > 1
                sum(length(parameter_mapping[q]) for q in 1:(b[1] - 1))
            else
                0
            end
            index = N + b[2]
            push!(remove, index)

            # need to now shuffle all the indices greater than this one down by 1
            for k in (b[2] + 1):length(parameter_mapping[b[1]])
                if parameter_mapping[b[1]][k] > parameter_number
                    parameter_mapping[b[1]][k] -= 1
                end
            end
            # and for subsequent models
            for j in (b[1] + 1):length(parameter_count)
                for k in 1:length(parameter_mapping[j])
                    if parameter_mapping[j][k] > parameter_number
                        parameter_mapping[j][k] -= 1
                    end
                end
            end
        end
    end

    sort!(unique!(remove))

    return parameter_mapping, remove
end

function _get_index_of_symbol(model::AbstractSpectralModel, symbol)::Int
    symbols = parameter_names(model)
    i = findfirst(==(symbol), symbols)
    if isnothing(i)
        error("Could not match symbol $symbol !")
    end
    return i
end

"""

Bind the symbols of `last(pair)` in all models indexes by `first(pair)`.
"""
function _bind_pairs!(prob::FittingProblem, pairs::Vararg{<:Pair{Int, Symbol}})
    binding = map(pairs) do pair
        model = prob.model.m[pair[1]]
        pair[1] => _get_index_of_symbol(model, pair[2])
    end |> collect
    _sort_binding!(binding)
    push!(prob.bindings, binding)
    return
end


"""
    bind!(prob::FittingProblem, pairs[, pairs...])

Bind parameters together within a [`FittingProblem`](@ref). Parameters bound
together will be mandated to have same value during the fit.

The binding must be specified in double or triple selector, which follow the
format:

```
pair := (model_index, :component_name, :parameter_symbol)
     := (model_index, :parameter_symbol)
```

The `model_index` is the index of the model in a multi-fit problem, i.e. `1`,
`2`, and so on.

The component name is a [`CompositeModel`](@ref) model name, e.g. `:a1`, or
`:c3`. This can be omitted if the model is not a [`CompositeModel`](@ref).

The parameter symbol is a symbol representing the field of the parameter in the
model. That is, `:K` or `:log10Flux`.

Bindings are specified using a chain of pairs `(root) => (target) [=>
(target)]`. The root parameter is kept as is, and all subsequent parameters are
bound to the root. Multiple chains of pairs may be specified in a single call
to `bind!`, or, alternatively, multiple bindings may be specified with
successive calls to `bind!`.

Bindings can be inspected with [`details`](@ref).

See also [`bindall!`](@ref).

## Examples

- Bind model 1's `K` parameter to model 2's second additive model's `K`:

  ```julia
  bind!(prob, (1, :K) => (2, :a2, :K))
  ```

- Bind model 3's `:a2.K` parameter to model4's `:m3.L` and model 6's `:a1.a`:

  ```julia
  bind!(prob, (3, :a2, :K) => (4, :m3, :K) => (6, :a1, :a))
  ```

Consider the following two models
```julia
model1 = PhotoelectricAbsorption() * (BlackBody() + PowerLaw())
model2 = PhotoelectricAbsorption() * (PowerLaw() + PowerLaw())

prob = FittingProblem(model1 => data1, model2 => data2)

# Bind the power law indices in the two models
bindall!(prob, :a)

# Bind the normalisation of powerlaws in the 2nd model:
bind!(prob, (2, :a1, :K) => (2, :a2, :K))

# To inspect the overall bindings.
details(prob)
```

!!! note
    Only free parameters can be bound together.
"""
function bind!(prob::FittingProblem, pairs::Vararg{<:Pair})
    for pair in pairs
        triples = _to_triple_vector(pair)
        bind!(prob, triples[1], triples[2:end])
    end
    return prob
end

function bind!(
        prob::FittingProblem,
        root::ParameterTriple,
        other::AbstractVector{<:ParameterTriple},
    )
    primary = parameter_index(prob, root)
    secondaries = map(i -> parameter_index(prob, i), other)
    existing = get(prob.bindings, primary, Int[])
    prob.bindings[primary] = unique!(vcat(existing, secondaries))
    return
end

"""
    bindall!(prob::FittingProblem, item[, item...])

Bind a common parameter across all models. The item is used to select the
parameter to bind, and may either be a single symbol, or a model-symbol double.

# Examples

```julia
# bind parameter `a` in all models
bindall!(prob, :a)

# bind parameter `K` in component `a3` in all models
bindall!(prob, (:a3, :K))

# multiple simultaneously
bindall!(prob, :E, (:a2, :K))
```
"""
function bindall!(
        prob::FittingProblem,
        items::Vararg{<:Union{<:Symbol, <:Tuple{Symbol, Symbol}}},
    )
    for item in items
        _item = if item isa Tuple
            item
        else
            (item,)
        end
        bind!(
            prob,
            ParameterTriple((1, _item...)),
            [ParameterTriple((i, _item...)) for i in 2:model_count(prob)],
        )
    end
    return prob
end

function _to_triple_vector(pair::Pair)
    out = [ParameterTriple(first(pair))]
    rem = last(pair)
    while (rem isa Pair)
        push!(out, ParameterTriple(first(rem)))
        rem = last(rem)
    end
    push!(out, ParameterTriple(rem))
    return out
end

parameter_index(prob::FittingProblem, t::ParameterTriple) = prob.lookup[t]


export bind!, bindall!, FittingProblem, FittableMultiModel, FittableMultiDataset, details
