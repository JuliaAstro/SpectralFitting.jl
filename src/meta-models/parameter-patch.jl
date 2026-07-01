"""
    ParameterPatchView{V<:AbstractVector}

Part of [`ParameterPatch`](@ref).

Used to create a view on a set of parameters that can be manipulated in an
intuitive way. Defines property setting / getting methods so that the parameters
associated with a particular model can be reached
```julia
p = ParameterPatchView(...)
p.a1 # returns a ModelPatchView
```

The possible properties that can be access are defined in a symbol vector.

Every "model" accessed this way will return a view on the parameters via
[`ModelPatchView`](@ref).
"""
struct ParameterPatchView{V <: AbstractVector}
    _values::V
    _symbols::Vector{Pair{Symbol, Symbol}}
    _access::Union{Nothing, BitVector}
    function ParameterPatchView(values::V, symbols; track = false) where {V <: AbstractVector}
        @assert size(values) == size(symbols)
        access = if track
            BitArray([false for i in symbols])
        else
            nothing
        end
        return new{V}(values, symbols, access)
    end
end

"""
    ModelPatchView{P<:ParameterPatchView}

Similar to [`ParameterPatchView`](@ref) but for accessing the parameter
components by symbol. Completes the ability to modify parameters as
```julia
pv = ParameterPatchView(...)
pv.a1.K = pv.a2.K * 3 + 2
```
"""
struct ModelPatchView{P <: ParameterPatchView}
    _parent::P
    _model::Symbol
end

Base.propertynames(pv::ParameterPatchView) = getfield(pv, :_symbols)
Base.propertynames(pv::ModelPatchView) = propertynames(getfield(pv, :_parent))

function Base.getproperty(pv::ModelPatchView, s::Symbol)
    parent = getfield(pv, :_parent)
    syms = getfield(parent, :_symbols)
    for (i, pair) in enumerate(syms)
        if first(pair) == getfield(pv, :_model) && last(pair) == s
            return getfield(parent, :_values)[i]
        end
    end
    error("Model $(getfield(pv, :_model)) has no parameter $(s)")
end

function Base.setproperty!(pv::ModelPatchView, s::Symbol, value)
    parent = getfield(pv, :_parent)
    syms = getfield(parent, :_symbols)
    for (i, pair) in enumerate(syms)
        if first(pair) == getfield(pv, :_model) && last(pair) == s
            if !isnothing(getfield(parent, :_access))
                getfield(parent, :_access)[i] = true
            end
            return getfield(parent, :_values)[i] = value
        end
    end
    error("Model $(getfield(pv, :_model)) has no parameter $(s)")
end

function Base.getproperty(pv::ParameterPatchView, s::Symbol)
    syms = getfield(pv, :_symbols)
    for pair in syms
        if first(pair) == s
            return ModelPatchView(pv, s)
        end
    end
    error("No such model component $(s)")
end

"""
    ParameterPatch(model; patch::Function)

An [`AbstractModelWrapper`](@ref) that can be used to manipulate the parameters
of the model it wraps. For example

```julia
model = PowerLaw() + PowerLaw()

function patcher!(p)
    # any arbitrary function may be defined here
    p.a1.K = 3 * p.a2.K + sqrt(p.a2.a)
end

patched_model = ParameterPatch(model; patch = patcher!)

# set the patched parameter as frozen so it is not fitted
# failing to do so may ruin a fit
patched_model.a1.K.frozen = true
```

When the model is invoked, it will call the `patch` function to manipulate the
parameters as desired.

!!! warning

    This wrapper is relatively new and not extensively tested. It should be
    noted nothing is done to validate if a parameter that you are modifying is
    actually a `frozen` parameter. This will change in the future, and
    parameters patched this way will be labelled as `bound`.

"""
struct ParameterPatch{M, T, K, F} <: AbstractModelWrapper{M, T, K}
    model::M
    symbols::Vector{Pair{Symbol, Symbol}}
    patch!::F
    function ParameterPatch(
            model::AbstractSpectralModel{T, K},
            symbols,
            patch!::F,
        ) where {T, K, F}
        return new{typeof(model), T, K, F}(model, symbols, patch!)
    end
end

function remake_with_parameters(
        model::ParameterPatch{T, K},
        params::AbstractVector,
    ) where {T, K}
    model.patch!(ParameterPatchView(params, model.symbols))
    return ParameterPatch(remake_with_parameters(model.model, params), model.symbols, model.patch!)
end

function ParameterPatch(model::AbstractSpectralModel{T, K}; patch::Function) where {T, K}
    _, syms = _all_parameters_with_symbols(model)
    model_copy = copy(model)
    # work out which parameters are being patched
    pp = ParameterPatch(model_copy, syms, patch)
    pvec, symbols = _all_parameters_with_symbols(model_copy)
    patch_view = ParameterPatchView(get_value.(pvec), symbols; track = true)
    pp.patch!(patch_view)
    for (i, patched) in enumerate(getfield(patch_view, :_access))
        if patched
            pvec[i].patched = true
        end
    end
    return pp
end

"""
    apply_patch!(model::ParameterPatch)

Apply a patch to the model (i.e. use the `patch!` function to update the model
parameters).
"""
function apply_patch!(model::ParameterPatch)
    pvec = parameter_vector(model.model)
    params = get_value.(pvec)
    model.patch!(ParameterPatchView(params, model.symbols))
    for (v, p) in zip(params, pvec)
        set_value!(p, v)
    end
    return model
end

function update_model!(model::ParameterPatch, result::FitResultSlice)
    update_model!(backing_model(model), result)
    return apply_patch!(model)
end

export ParameterPatch, ParameterPatchView, apply_patch!
