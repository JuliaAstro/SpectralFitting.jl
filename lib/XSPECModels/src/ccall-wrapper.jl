UNTRACKED_ERROR_F64 = zeros(Float64, 1)
UNTRACKED_ERROR_F32 = zeros(Float32, 1)

DEFAULT_MODEL_FLOAT_TYPE = Float64
DEFAULT_MODEL_INT_TYPE = Cint

@inline function get_untracked_error(T::Type)
    if T === Float64
        UNTRACKED_ERROR_F64
    elseif T === Float32
        UNTRACKED_ERROR_F32
    else
        error("Unsupported type $(T)")
    end
end

resize_untracked_error!(_error, n) = resize!(_error, n)

_unsafe_unwrap_parameter(ptr, N, ::AbstractSpectralModelKind, T::Type) =
    unsafe_wrap(Vector{T}, ptr, N, own = false)
# additive models don't include normalisation, so we offset
_unsafe_unwrap_parameter(ptr, N, ::Additive, T::Type) =
    unsafe_wrap(Vector{T}, ptr + sizeof(T), N - 1, own = false)

@inline function unsafe_parameter_vector(
    model_ref::Ref{M},
) where {M<:AbstractSpectralModel{T}} where {T}
    N = length(fieldnames(M))
    ptr = convert(Ptr{T}, pointer_from_objref(model_ref))
    # reinterpret as an unowned vector of `T`
    _unsafe_unwrap_parameter(ptr, N, modelkind(M), T)
end

# pointer hackery
function unsafe_parameter_vector_conditioned(
    model_ref::Ref{<:AbstractSpectralModel{T}},
) where {T}
    if implementation(typeof(model_ref[])) == JuliaImplementation()
        throw("This method is only for `XSPECImplementation` models.")
    end
    if !isbitstype(T)
        throw("Can only reinterpret to bits type (convert `$T`).")
    end
    unsafe_parameter_vector(model_ref)
end

"""
    function _unsafe_ffi_invoke!(
        output,
        error_vec,
        input,
        params,
        ModelType::Type{<:AbstractSpectralModel},
    )

Wrapper to do a foreign function call to an XSPEC model.

See also [`_safe_ffi_invoke!`](@ref).
"""
function _unsafe_ffi_invoke!(
    output,
    error_vec,
    input,
    params,
    ModelType::Type{<:AbstractSpectralModel},
)
    error("""
          Not yet implemented for $(ModelType).
          If you set `generate_ffi_call=false` when using `@xspecmodel`,
          you must implement this function. Consult the documentation of
          `@xspecmodel` for reference.
          """)
end

"""
    function _safe_ffi_invoke!(output, input, params, ModelType::Type{<:AbstractSpectralModel})

Wrapper to do a foreign function call to an XSPEC model.
"""
function _safe_ffi_invoke!(output, input, params, ModelType::Type{<:AbstractSpectralModel})
    error(
        "This error should be unreachable. Please open an issue with your use-case and include the stack trace.",
    )
end

macro wrap_xspec_model_ccall(
    FloatType,
    IntType,
    func_name,
    callsite,
    input,
    output,
    err,
    parameters,
    spec_number,
    init_string,
)
    quote
        ccall(
            ($(func_name), $(callsite)),
            Cvoid,
            (
                Ptr{$(FloatType)},
                $(IntType),
                Ptr{$(FloatType)},
                $(IntType),
                Ptr{$(FloatType)},
                Ptr{$(FloatType)},
                Cstring,
            ),
            $(input),
            length($(input)) - 1,
            $(parameters),
            $(spec_number),
            $(output),
            $(err),
            $(init_string),
        )
    end |> esc
end

"""
    @xspecmodel [type=Float64] [ff_call_site] model

Used to wrap additional XSPEC models, generating the needed [`AbstractSpectralModel`](@ref)
implementation.

The `type` keyword specifies the underlying type to coerce input and output arrays to, as different
implementations may have incompatible number of bits. The `ff_call_site` is the foreign function call site,
which is the first argument to `ccall`, and follows the same conventions. The `model` is a struct, which must subtype
[`AbstractSpectralModel`](@ref).

If the callsite is not specified, the user must implement [`_unsafe_ffi_invoke!`](@ref).

# Examples

```julia
@xspecmodel :C_powerlaw struct XS_PowerLaw{T} <: AbstractSpectralModel{T, Additive}
    "Normalisation."
    K::T
    "Photon index."
    a::T
end

# constructor has default values
function XS_PowerLaw(; K = FitParam(1.0), a = FitParam(1.0))
    XS_PowerLaw{typeof(K)}(K, a)
end
```

We define a new structure `XS_PowerLaw` with two parameters, but since the model is [`Additive`](@ref),
only a single parameter (`a`) is passed to the XSPEC function. The function we bind to this model
is `:C_powerlaw` from the XSPEC C wrappers.

The macro will then generate the following functions
- [`implementation`](@ref)
- [`invoke!`](@ref)
- [`_safe_ffi_invoke!`](@ref)

If a callsite was specified, it will also generate:
- [`_unsafe_ffi_invoke!`](@ref)
"""
macro xspecmodel(args...)
    # parse args
    model = args[end]
    c_function = nothing
    model_float_type::Type = DEFAULT_MODEL_FLOAT_TYPE
    model_int_type::Type = DEFAULT_MODEL_INT_TYPE

    if (length(args) > 1)
        for arg in args[1:(end-1)]
            if (arg isa QuoteNode) || (arg.head !== :(=))
                # parse only optional positional
                if isnothing(c_function)
                    c_function = arg
                else
                    error("Too many positionals.")
                end
                continue
            end
            # parse key value pair
            key = first(arg.args)
            value = last(arg.args)
            if key == :type
                model_float_type = if value == :Float32
                    Float32
                elseif value == :Float64
                    Float64
                else
                    error("Unsupported type $(value)")
                end
            else
                error("Unknown key value pair: $(arg)")
            end
        end
    end

    model_args = model.args[3]
    model_name = model.args[2].args[1].args[1]
    model_type_params = model.args[2].args[1].args[2:end]
    model_kind = model.args[2].args[end]
    symbols = [i.args[1] for i in model_args.args if (i isa Expr) && (i.head == :(::))]

    if model_kind.args[2] == :Additive
        # get rid of the normalisation from function arguments
        symbols = symbols[2:end]
    end

    _ffi_type_guard = _build_ffi_type_guard(model_name, model_float_type)

    _unsafe_call_def = if !isnothing(c_function)
        if c_function isa QuoteNode
            call_symbol = c_function
            call_site = libXSFunctions
        else
            call_symbol = c_function.args[1]
            call_site = c_function.args[2]
        end
        _build_unsafe_ffi_call(
            model_name,
            call_symbol,
            call_site,
            model_float_type,
            model_int_type,
        )
    else
        :()
    end

    quote
        Base.@__doc__ struct $(model_name){$(model_type_params...)} <: $(model_kind)
            $(model_args.args...)
        end

        SpectralFitting.implementation(::Type{<:$(model_name)}) = XSPECImplementation()

        function SpectralFitting.invoke!(output, input, model::$(model_name))
            # allocate the model
            model_ref = Ref(model)
            params = @noinline XSPECModels.unsafe_parameter_vector_conditioned(model_ref)
            XSPECModels._safe_ffi_invoke!(vec(output), input, params, typeof(model))
        end

        $(_ffi_type_guard)

        @inline function XSPECModels._safe_ffi_invoke!(
            output::AbstractVector{<:$(model_float_type)},
            input::AbstractVector{<:$(model_float_type)},
            params::AbstractVector{<:$(model_float_type)},
            ModelType::Type{<:$(model_name)},
        )
            @assert length(output) + 1 == length(input)
            SpectralFitting.ensure_model_data($(model_name))

            error_vec = XSPECModels.get_untracked_error($(model_float_type))
            if length(error_vec) < length(output)
                XSPECModels.resize_untracked_error!(error_vec, length(output))
            end
            XSPECModels._unsafe_ffi_invoke!(output, error_vec, input, params, ModelType)
        end

        $(_unsafe_call_def)

        $(model_name)
    end |> esc
end

function _build_ffi_type_guard(model_name, model_float_type)
    if model_float_type != DEFAULT_MODEL_FLOAT_TYPE
        :(@inline function XSPECModels._safe_ffi_invoke!(
            output,
            input,
            params,
            ModelType::Type{<:$(model_name)},
        )
            f_typed_output = convert.($(model_float_type), (output))
            f_typed_input = convert.($(model_float_type), (input))
            f_typed_params = convert.($(model_float_type), (params))
            XSPECModels._safe_ffi_invoke!(
                f_typed_output,
                f_typed_input,
                f_typed_params,
                ModelType,
            )
            @. output = f_typed_output
        end)
    else
        :()
    end
end

function _build_unsafe_ffi_call(
    model_name,
    call_symbol,
    call_site,
    model_float_type,
    model_int_type,
)
    :(
        function XSPECModels._unsafe_ffi_invoke!(
            output,
            error_vec,
            input,
            params,
            ::Type{<:$(model_name)},
        )
            XSPECModels.@wrap_xspec_model_ccall(
                $(model_float_type),
                $(model_int_type),
                $(call_symbol),
                $(call_site),
                input,
                output,
                error_vec,
                params,
                1,
                ""
            )
        end
    )
end


export @xspecmodel, @wrap_xspec_model_ccall
