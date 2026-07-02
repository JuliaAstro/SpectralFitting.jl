export AbstractCompositeOperator,
    AdditionOperator,
    MultiplicationOperator,
    ConvolutionOperator,
    operation_symbol,
    CompositeModel

"""
    abstract type AbstractCompositeOperator

Superype of all composition operators. Used to implement the model algebra of
[`AbstractSpectralModelKind`](@ref) through a trait system.

The Julia symbol corresponding to a given `AbstractCompositeOperator` may be
obtained through [`operation_symbol`](@ref).
"""
abstract type AbstractCompositeOperator end
"""
    AdditionOperator <: AbstractCompositeOperator
    AdditionOperator()

Corresponds to the `:(+)` symbol.
"""
struct AdditionOperator <: AbstractCompositeOperator end
"""
    MultiplicationOperator <: AbstractCompositeOperator
    MultiplicationOperator()

Corresponds to the `:(*)` symbol.
"""
struct MultiplicationOperator <: AbstractCompositeOperator end
"""
    ConvolutionOperator <: AbstractCompositeOperator
    ConvolutionOperator()

Has no corresponding symbol, since it invokes a function call `C(A)`.
"""
struct ConvolutionOperator <: AbstractCompositeOperator end

"""
    operation_symbol(::AbstractCompositeOperator)
    operation_symbol(::Type{<:AbstractCompositeOperator})

Obtain the model symbol from a given [`AbstractCompositeOperator`](@ref).
"""
operation_symbol(::Type{<:AdditionOperator}) = :(+)
operation_symbol(::Type{<:MultiplicationOperator}) = :(*)
operation_symbol(::O) where {O <: AbstractCompositeOperator} = operation_symbol(O)

combine_components!(::Type{<:AdditionOperator}, left, right) = @. right = left + right
combine_components!(::Type{<:MultiplicationOperator}, left, right) = @. right = left * right

# algebra grammar
add_models(_, _, ::M1, ::M2) where {M1, M2} =
    throw("Left and right models must be Additive.")
add_models(m1, m2, ::Additive, ::Additive) = CompositeModel(m1, m2, AdditionOperator())
add_models(m1::M1, m2::M2) where {M1, M2} = add_models(m1, m2, modelkind(M1), modelkind(M2))
Base.:+(m1::M1, m2::M2) where {M1 <: AbstractSpectralModel, M2 <: AbstractSpectralModel} =
    add_models(m1, m2, modelkind(M1), modelkind(M2))

mult_models(_, _, ::M1, ::M2) where {M1, M2} = throw("Left model must be Multiplicative.")
mult_models(m1, m2, ::Multiplicative, ::AbstractSpectralModelKind) =
    CompositeModel(m1, m2, MultiplicationOperator())
mult_models(m1::M1, m2::M2) where {M1, M2} =
    mult_models(m1, m2, modelkind(M1), modelkind(M2))
Base.:*(m1::M1, m2::M2) where {M1 <: AbstractSpectralModel, M2 <: AbstractSpectralModel} =
    mult_models(m1, m2, modelkind(M1), modelkind(M2))

conv_models(_, _, ::M1, ::M2) where {M1, M2} =
    throw("Left model must be Convolutional and right model must be Additive.")
conv_models(m1, m2, ::Convolutional, ::Additive) =
    CompositeModel(m1, m2, ConvolutionOperator())
conv_models(m1::M1, m2::M2) where {M1, M2} =
    conv_models(m1, m2, modelkind(M1), modelkind(M2))
(m1::AbstractSpectralModel)(m2::M2) where {M2 <: AbstractSpectralModel} = conv_models(m1, m2)

"""
    CompositeModel{T,K,Op,M1,M2} <: AbstractSpectralModel{T,K}
    CompositeModel(left_model, right_model, op::AbstractCompositeOperator)

Type resulting from operations combining any number of
[`AbstractSpectralModel`](@ref) via the model algebra defined from
[`AbstractSpectralModelKind`](@ref).

```julia
model = PhotoelectricAbsorption() * (PowerLaw() + BlackBody())
typeof(model) <: CompositeModel # true
```

Each operation binary operation in the model algebra is encoded in the
parametric types of the `CompositeModel`, where the operation is given by an
[`AbstractCompositeOperator`](@ref).  Composite models adopt the model kind of
the `right` model, i.e. `M2`, and obey the model algebra accordingly.

Composite models very rarely need to be constructed directly, and are instead
obtained by regular model operations.

The `propertynames` and `getproperty` methods for `CompositeModel` are
overwritten to access all models in the model tree. For example:

```julia
julia> propertynames(PowerLaw() + DeltaLine())
(:a1, :a2)
```

`CompositeModel` has [`destructure`](@ref) for working abstractly with model
trees.
"""
struct CompositeModel{T, K, Operator, M1 <: AbstractSpectralModel, M2 <: AbstractSpectralModel} <:
    AbstractSpectralModel{T, K}
    left::M1
    right::M2
end

function CompositeModel(
        left::AbstractSpectralModel{T},
        right::AbstractSpectralModel{T, K},
        ::Operator;
        make_copies = true,
    ) where {T, K, Operator}
    if make_copies
        return CompositeModel{T, K, Operator, typeof(left), typeof(right)}(copy(left), copy(right))
    else
        return CompositeModel{T, K, Operator, typeof(left), typeof(right)}(left, right)
    end
end

function Base.copy(m::CompositeModel{T, K, Op}) where {T, K, Op}
    return CompositeModel(getfield(m, :left), getfield(m, :right), Op())
end

is_composite(::M) where {M <: AbstractSpectralModel} = is_composite(M)
is_composite(::Type{<:CompositeModel}) = true
is_composite(::Type{<:AbstractSpectralModel}) = false

function implementation(::Type{<:CompositeModel{T, K, Op, M1, M2}}) where {T, K, Op, M1, M2}
    if implementation(M1) isa XSPECImplementation || implementation(M2) isa XSPECImplementation
        return XSPECImplementation()
    else
        return JuliaImplementation()
    end
end

objective_cache_count(model::AbstractSpectralModel) = 1
objective_cache_count(model::CompositeModel) =
    objective_cache_count(getfield(model, :left)) +
    objective_cache_count(getfield(model, :right))

parameter_vector(model::CompositeModel) = vcat(
    parameter_vector(getfield(model, :left)),
    parameter_vector(getfield(model, :right)),
)

parameter_count(model::CompositeModel) =
    parameter_count(getfield(model, :left)) + parameter_count(getfield(model, :right))

function remake_with_parameters(model::CompositeModel{T, K, Op}, params::Tuple) where {T, K, Op}
    P_left = parameter_count(getfield(model, :left))
    return CompositeModel(
        remake_with_parameters(getfield(model, :left), params[1:P_left]),
        remake_with_parameters(getfield(model, :right), params[(P_left + 1):end]),
        Op();
        make_copies = false,
    )
end

function _inner_invokemodel!(
        outputs,
        domain,
        model::CompositeModel{<:Number, K, Op, M1, M2},
    ) where {K, Op, M1, M2}
    @assert size(outputs, 2) >= objective_cache_count(model) "Too few flux arrays allocated for this model ($(size(outputs)))."
    RetType = typeof(@views(outputs[:, 1]))

    return @views begin
        _out_right = if is_composite(M2)
            outputs[:, 1:end]
        else
            outputs[:, 1:1]
        end

        _inner_invokemodel!(_out_right, domain, getfield(model, :right))
        if Op <: ConvolutionOperator
            # Use the same principle output as the right hand side
            _out_left = if is_composite(M1)
                outputs[:, 1:end]
            else
                outputs[:, 1:1]
            end
            return _inner_invokemodel!(_out_left, domain, getfield(model, :left))
        else
            _out_left = if is_composite(M1)
                outputs[:, 2:end]
            else
                outputs[:, 2:2]
            end
            _inner_invokemodel!(_out_left, domain, getfield(model, :left))
            return combine_components!(Op, outputs[:, 2], outputs[:, 1])
        end
    end
end

# printing
struct ModelKindCounts{N}
    additive::Int
    multiplicative::Int
    convolutional::Int
    syms::NTuple{N, Symbol}
end

function destructure(::Type{<:AbstractSpectralModel{T, K}}, d::ModelKindCounts) where {T, K}
    if K == Additive
        s = Symbol("a$(d.additive + 1)")
        m = ModelKindCounts(d.additive + 1, d.multiplicative, d.convolutional, (d.syms..., s))
    elseif K == Multiplicative
        s = Symbol("m$(d.multiplicative + 1)")
        m = ModelKindCounts(d.additive, d.multiplicative + 1, d.convolutional, (d.syms..., s))
    elseif K == Convolutional
        s = Symbol("c$(d.convolutional + 1)")
        m = ModelKindCounts(d.additive, d.multiplicative, d.convolutional + 1, (d.syms..., s))
    else
        error("Unknown model kind: $K")
    end
    return (s, m)
end

function destructure(
        ::Type{<:CompositeModel{T, K, Op, Left, Right}},
        d::ModelKindCounts,
    ) where {T, K, Op, Left, Right}
    exp_left, d_left = destructure(Left, d)
    exp_right, d_right = destructure(Right, d_left)
    if Op <: ConvolutionOperator
        return :($(exp_left)($exp_right)), d_right
    else
        op = SpectralFitting.operation_symbol(Op)
        return :($(op)($exp_left, $exp_right)), d_right
    end
end

destructure(M::Type{<:CompositeModel}) = destructure(M, ModelKindCounts(0, 0, 0, ()))
flatten_models(m::AbstractSpectralModel) = (m,)

function flatten_models(m::CompositeModel)
    right = flatten_models(getfield(m, :right))
    left = flatten_models(getfield(m, :left))
    return (left..., right...)
end

struct DestructuredCompositeModel{M}
    expr::String
    models::Vector{Pair{Symbol, M}}
end

function destructure(@nospecialize(m::CompositeModel))
    expr, info = destructure(typeof(m))
    all_models = flatten_models(m)
    return DestructuredCompositeModel(
        "$(expr)",
        Pair{Symbol, AbstractSpectralModel}[k => v for (k, v) in zip(info.syms, all_models)],
    )
end

function _printinfo(
        io::IO,
        @nospecialize(model::CompositeModel);
        bindings = nothing,
        name = "CompositeModel",
    )
    expr_buffer = 3
    sym_buffer = 6

    dis = destructure(model)

    println(io, "$name with $(length(dis.models)) model components:")
    print(io, " "^expr_buffer)
    printstyled(io, dis.expr, color = :cyan)
    println(io)
    println(io, "Model key and parameters:")

    param_index = 1
    for (sym, m) in dis.models
        basename = Base.typename(typeof(m)).name

        printstyled(io, lpad("$sym", sym_buffer), color = :cyan)
        print(io, " => ")

        buff = IOBuffer()
        _printinfo(
            IOContext(buff, io),
            m;
            parameter_indent = 3,
            bindings = if isnothing(bindings)
                nothing
            else
                get(bindings, sym, nothing)
            end,
        )
        s = String(take!(buff))
        println(io, strip(indent(s, 4)))
    end
    return
end

# explicitly disable interface for ConstructionBase.jl
ConstructionBase.setproperties(::CompositeModel, ::NamedTuple) =
    throw("Cannot be used with `CompositeModel`.")
ConstructionBase.constructorof(::Type{<:CompositeModel}) =
    throw("Cannot be used with `CompositeModel`.")

function Base.propertynames(model::CompositeModel)
    _, info = destructure(typeof(model))
    return info.syms
end

function Base.getproperty(model::CompositeModel, symb::Symbol)
    model_map = destructure(model)
    ind = findfirst(i -> i[1] == symb, model_map.models)
    if isnothing(ind)
        error("$(_model_name(model)) has no component $symb")
    else
        return last(model_map.models[ind])
    end
end

function parameter_names(::Type{<:CompositeModel{T, K, O, M1, M2}}) where {T, K, O, M1, M2}
    return (parameter_names(M1)..., parameter_names(M2)...)
end

function _all_parameters_with_symbols(model::CompositeModel{T}) where {T}
    dis = destructure(model)
    all_params = T[]
    all_syms = Pair{Symbol, Symbol}[]
    for (model_sym, model) in dis.models
        params, syms = _all_parameters_with_symbols(model)
        append!(all_params, params)
        append!(all_syms, [model_sym => i[2] for i in syms])
    end
    return all_params, all_syms
end
