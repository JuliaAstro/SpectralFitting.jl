# =============================================================================
# xspec_string.jl
#
# Bidirectional translation between XSPEC model strings and XSPECModels.jl
# composite model expressions.
#
# Issue: https://github.com/JuliaAstro/SpectralFitting.jl/issues/187
#
# PUBLIC API
# ──────────────────────────────────────────────────────────────────────────────
#   xspec"phabs*powerlaw"              → instantiated CompositeModel (macro)
#   parse_xspec_model_string(str)      → Julia Expr (runtime, no eval needed)
#   xspec_model_string(model)          → XSPEC string from any model
#   XSPEC_MODEL_NAMES                  → Dict{String,Symbol} of all mappings
#
# ONLY MODELS THAT ACTUALLY EXIST IN THIS PACKAGE ARE LISTED.
# Verified against the live source of:
#   lib/XSPECModels/src/additive.jl
#   lib/XSPECModels/src/multiplicative.jl
#   lib/XSPECModels/src/convolutional.jl
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# 1.  FORWARD LOOKUP TABLE
#     XSPEC model-string token  →  Julia type Symbol
#
#     Every entry has been cross-checked against @xspecmodel declarations in
#     the source files.  Only types that are actually exported by this package
#     appear here.
# ─────────────────────────────────────────────────────────────────────────────

"""
    XSPEC_MODEL_NAMES :: Dict{String, Symbol}

Maps every recognised XSPEC model-string token to the corresponding
XSPECModels.jl Julia type name.

Only models that are actually defined and exported by this package are
included.  Both canonical XSPEC names and common abbreviations ("po",
"gaus", "bb") are supported; all keys are lowercase.

```julia
XSPEC_MODEL_NAMES["powerlaw"]   # :XS_PowerLaw
XSPEC_MODEL_NAMES["tbabs"]      # :XS_NeutralHydrogenAbsorption
XSPEC_MODEL_NAMES["cflux"]      # :XS_CalculateFlux
```
"""
const XSPEC_MODEL_NAMES = Dict{String,Symbol}(

    # ── Additive models (from additive.jl) ───────────────────────────────────
    # @xspecmodel :C_powerlaw → XS_PowerLaw
    "powerlaw"   => :XS_PowerLaw,
    "po"         => :XS_PowerLaw,          # common XSPEC abbreviation

    # @xspecmodel :C_zcutoffpl → XS_CutOffPowerLaw
    "cutoffpl"   => :XS_CutOffPowerLaw,
    # "zcutoffpl" is the C function prefix; canonical XSPEC name is "cutoffpl"
    # "zcutoffpl"  => :XS_CutOffPowerLaw,

    # @xspecmodel :C_bbody → XS_BlackBody
    "bbody"      => :XS_BlackBody,
    "bb"         => :XS_BlackBody,         # common XSPEC abbreviation

    # @xspecmodel :C_bremss → XS_BremsStrahlung
    "bremss"     => :XS_BremsStrahlung,

    # @xspecmodel :C_kerrdisk → XS_KerrDisk
    "kerrdisk"   => :XS_KerrDisk,

    # @xspecmodel :C_kyrline → XS_KyrLine
    "kyrline"    => :XS_KyrLine,

    # @xspecmodel :C_laor → XS_Laor
    "laor"       => :XS_Laor,

    # @xspecmodel :C_diskline → XS_DiskLine
    "diskline"   => :XS_DiskLine,

    # @xspecmodel :C_gaussian → XS_Gaussian
    "gaussian"   => :XS_Gaussian,
    "gaus"       => :XS_Gaussian,          # common XSPEC abbreviation

    # @xspecmodel :C_jet → XS_Jet
    "jet"        => :XS_Jet,

    # @xspecmodel :C_optxagnf → XS_Optxagnf
    "optxagnf"   => :XS_Optxagnf,

    # @xspecmodel :C_diskbb → XS_DiskBlackBody
    "diskbb"     => :XS_DiskBlackBody,

    # ── Multiplicative models (from multiplicative.jl) ────────────────────────
    # @xspecmodel :C_phabs → XS_PhotoelectricAbsorption
    "phabs"      => :XS_PhotoelectricAbsorption,

    # @xspecmodel :C_wndabs → XS_WarmAbsorption
    "wndabs"     => :XS_WarmAbsorption,

    # @xspecmodel :C_tbabs → XS_NeutralHydrogenAbsorption
    "tbabs"      => :XS_NeutralHydrogenAbsorption,

    # ── Convolutional models (from convolutional.jl) ──────────────────────────
    # @xspecmodel :C_cflux → XS_CalculateFlux
    "cflux"      => :XS_CalculateFlux,

    # @xspecmodel :C_kerrconv → XS_Kerrconv
    "kerrconv"   => :XS_Kerrconv,
)


# ─────────────────────────────────────────────────────────────────────────────
# 2.  REVERSE LOOKUP TABLE
#     Julia type Symbol  →  canonical XSPEC token
#
#     Built automatically from XSPEC_MODEL_NAMES.
#     Where multiple tokens map to the same type (aliases like "po"/"powerlaw"),
#     the LONGEST token is chosen as the canonical form because it is more
#     human-readable in reconstructed strings.
# ─────────────────────────────────────────────────────────────────────────────

const _XSPEC_REVERSE = let
    d = Dict{Symbol,String}()
    for (xname, jtype) in XSPEC_MODEL_NAMES
        # keep the longest XSPEC name as canonical
        if !haskey(d, jtype) || length(xname) > length(d[jtype])
            d[jtype] = xname
        end
    end
    d
end


# ─────────────────────────────────────────────────────────────────────────────
# 3.  TOKENISER
# ─────────────────────────────────────────────────────────────────────────────

@enum _TokKind::UInt8 begin
    _TOK_NAME      # "powerlaw", "phabs", etc.
    _TOK_STAR      # *
    _TOK_PLUS      # +
    _TOK_LPAREN    # (
    _TOK_RPAREN    # )
    _TOK_EOF
end

struct _Token
    kind  :: _TokKind
    value :: String   # non-empty only for _TOK_NAME
end

"""
Internal: tokenise an XSPEC model string into `_Token` objects.
Lowercases the whole input so the parser is case-insensitive.
"""
function _tokenize(src::AbstractString)
    tokens = _Token[]
    s = strip(lowercase(src))
    i = firstindex(s)
    while i <= lastindex(s)
        c = s[i]
        if isspace(c)
            i = nextind(s, i)
        elseif c == '*'
            push!(tokens, _Token(_TOK_STAR,   "")); i = nextind(s, i)
        elseif c == '+'
            push!(tokens, _Token(_TOK_PLUS,   "")); i = nextind(s, i)
        elseif c == '('
            push!(tokens, _Token(_TOK_LPAREN, "")); i = nextind(s, i)
        elseif c == ')'
            push!(tokens, _Token(_TOK_RPAREN, "")); i = nextind(s, i)
        elseif isletter(c) || isdigit(c)
            j = i
            while j <= lastindex(s) && (isletter(s[j]) || isdigit(s[j]) || s[j] == '_')
                j = nextind(s, j)
            end
            push!(tokens, _Token(_TOK_NAME, s[i:prevind(s, j)]))
            i = j
        else
            error("XSPEC string parse error: unexpected character '$c' in \"$src\"")
        end
    end
    push!(tokens, _Token(_TOK_EOF, ""))
    return tokens
end


# ─────────────────────────────────────────────────────────────────────────────
# 4.  RECURSIVE-DESCENT PARSER  →  Julia Expr
#
#   Grammar (operator precedence matches XSPEC: * binds tighter than +):
#
#       expr   ::= term  ( '+' term  )*
#       term   ::= factor( '*' factor)*
#       factor ::= '(' expr ')' | NAME
#
#   This mirrors how XSPEC itself parses its model strings and is the correct
#   astrophysical interpretation:
#       "phabs*powerlaw+bbody"   →  (phabs*powerlaw) + bbody
#       "phabs*(powerlaw+bbody)" →  phabs * (powerlaw+bbody)
# ─────────────────────────────────────────────────────────────────────────────

mutable struct _Parser
    tokens :: Vector{_Token}
    pos    :: Int
    src    :: String    # kept only for error messages
end

@inline _peek(p::_Parser)    = p.tokens[p.pos]
@inline function _advance!(p::_Parser)
    t = p.tokens[p.pos]; p.pos += 1; t
end
function _expect!(p::_Parser, k::_TokKind)
    t = _advance!(p)
    t.kind == k && return t
    error(
        "XSPEC string parse error: expected $k but got $(t.kind) " *
        "('$(t.value)') in \"$(p.src)\""
    )
end

function _parse_expr!(p::_Parser)
    left = _parse_term!(p)
    while _peek(p).kind == _TOK_PLUS
        _advance!(p)
        right = _parse_term!(p)
        left  = Expr(:call, :+, left, right)
    end
    return left
end

function _parse_term!(p::_Parser)
    left = _parse_factor!(p)
    while _peek(p).kind == _TOK_STAR
        _advance!(p)
        right = _parse_factor!(p)
        left  = Expr(:call, :*, left, right)
    end
    return left
end

function _parse_factor!(p::_Parser)
    t = _peek(p)

    if t.kind == _TOK_LPAREN
        _advance!(p)
        inner = _parse_expr!(p)
        _expect!(p, _TOK_RPAREN)
        return inner

    elseif t.kind == _TOK_NAME
        _advance!(p)
        name = t.value
        if haskey(XSPEC_MODEL_NAMES, name)
            return Expr(:call, XSPEC_MODEL_NAMES[name])   # e.g. :(XS_PowerLaw())
        else
            # Helpful prefix hint
            prefix  = name[1:min(3, ncodeunits(name))]
            similar = sort!(filter(k -> startswith(k, prefix), collect(keys(XSPEC_MODEL_NAMES))))
            hint    = isempty(similar) ? "" :
                      "\n  Did you mean one of: $(join(similar[1:min(5,end)], ", "))?"
            error(
                "XSPEC string parse error: unknown model name '$name' in " *
                "\"$(p.src)\".$hint\n" *
                "  See `XSPEC_MODEL_NAMES` for the full list of supported models.\n" *
                "  Note: only models implemented in XSPECModels.jl are supported."
            )
        end

    elseif t.kind == _TOK_EOF
        error("XSPEC string parse error: unexpected end of string in \"$(p.src)\"")

    else
        error(
            "XSPEC string parse error: unexpected token $(t.kind) " *
            "('$(t.value)') in \"$(p.src)\""
        )
    end
end


# ─────────────────────────────────────────────────────────────────────────────
# 5.  PUBLIC API — FORWARD  (XSPEC string → Julia)
# ─────────────────────────────────────────────────────────────────────────────

"""
    parse_xspec_model_string(s::AbstractString) → Expr

Parse an XSPEC model string and return a Julia `Expr` that, when `eval`'d,
constructs the equivalent composite model with default parameters.

This is the **runtime** version.  For a compile-time macro with zero overhead
use the `xspec"..."` string macro instead.

## Grammar

| Token        | Meaning                                  |
|:-------------|:-----------------------------------------|
| `name`       | Instantiate the named XSPEC model        |
| `m1 * m2`    | Multiplicative: `m1` acting on `m2`     |
| `a1 + a2`    | Additive: sum of two additive components |
| `( … )`      | Grouping (override default precedence)   |

`*` binds more tightly than `+`, matching XSPEC conventions.
The string is case-insensitive and whitespace-tolerant.

## Examples

```julia
parse_xspec_model_string("powerlaw")
# :(XS_PowerLaw())

parse_xspec_model_string("phabs*powerlaw")
# :(XS_PhotoelectricAbsorption() * XS_PowerLaw())

parse_xspec_model_string("tbabs*(powerlaw+bbody)")
# :(XS_NeutralHydrogenAbsorption() * (XS_PowerLaw() + XS_BlackBody()))

# Precedence: * binds tighter than +
parse_xspec_model_string("phabs*powerlaw+bbody")
# :((XS_PhotoelectricAbsorption() * XS_PowerLaw()) + XS_BlackBody())
```

See also [`XSPEC_MODEL_NAMES`](@ref) for the full list of supported names and
[`xspec_model_string`](@ref) for the reverse direction.
"""
function parse_xspec_model_string(s::AbstractString)
    isempty(strip(s)) && error("XSPEC string parse error: empty string provided")
    tokens = _tokenize(s)
    p      = _Parser(tokens, 1, String(s))
    expr   = _parse_expr!(p)
    if _peek(p).kind != _TOK_EOF
        leftover = _peek(p).value
        error(
            "XSPEC string parse error: unexpected token '$leftover' at position $(p.pos)" *
            " in \"$s\"\n  Check for unmatched parentheses or double operators."
        )
    end
    return expr
end


"""
    @xspec_str(s) → model

String macro: parse an XSPEC model string **at compile time** and return an
instantiated SpectralFitting.jl composite model with default parameters.

The macro expands during compilation, so there is **zero runtime overhead**
beyond constructing the model structs themselves.

```julia
using XSPECModels

# Simple absorbed power-law (most common AGN baseline)
m = xspec"phabs*powerlaw"

# Thermal plasma with neutral hydrogen absorption
m = xspec"tbabs*powerlaw"

# Two-component additive model with absorption
m = xspec"phabs*(powerlaw+bbody)"

# Multi-component with convolution
m = xspec"cflux*powerlaw"

# Precedence: * tighter than +
#   "phabs*powerlaw+bbody"  →  (phabs*powerlaw) + bbody
m = xspec"phabs*powerlaw+bbody"
```

The string is case-insensitive and whitespace-tolerant.

!!! note "Supported models"
    Only models actually implemented in XSPECModels.jl are supported.
    See [`XSPEC_MODEL_NAMES`](@ref) for the full list.

See also [`parse_xspec_model_string`](@ref) for the runtime version and
[`xspec_model_string`](@ref) for the reverse direction.
"""
macro xspec_str(s)
    return esc(parse_xspec_model_string(s))
end


# ─────────────────────────────────────────────────────────────────────────────
# 6.  PUBLIC API — REVERSE  (model → XSPEC string)
#
#     Works for both atomic XSPECModels types and SpectralFitting.jl
#     CompositeModel trees built from them.
#
#     CompositeModel{T,K,Operator,M1,M2} has:
#       - field :left  → M1 instance
#       - field :right → M2 instance
#       - type param 3 (Operator) encodes the operation:
#           AdditionOperator       → "+"
#           MultiplicationOperator → "*"
#           ConvolutionOperator    → treated as "*" for string purposes
# ─────────────────────────────────────────────────────────────────────────────

"""
    xspec_model_string(model) → String

Convert a SpectralFitting.jl model back to its XSPEC model string.
Supports atomic XSPECModels types and any `CompositeModel` built from them.

```julia
xspec_model_string(XS_PowerLaw())
# "powerlaw"

xspec_model_string(XS_PhotoelectricAbsorption() * XS_PowerLaw())
# "phabs*powerlaw"

xspec_model_string(XS_NeutralHydrogenAbsorption() * (XS_PowerLaw() + XS_BlackBody()))
# "tbabs*(powerlaw+bbody)"

xspec_model_string(
    XS_PhotoelectricAbsorption() * (XS_Gaussian() * (XS_PowerLaw() + XS_BlackBody()))
)
# "phabs*(gaussian*(powerlaw+bbody))"
```

Throws `ArgumentError` if the model contains a type not in [`XSPEC_MODEL_NAMES`](@ref).

See also [`parse_xspec_model_string`](@ref), [`xspec"..."`](@ref @xspec_str).
"""
xspec_model_string(model) = _to_xspec(model, false)


# ── Dispatch: atomic model ────────────────────────────────────────────────────
function _to_xspec(model::AbstractSpectralModel, needs_parens::Bool)
    T     = typeof(model)
    tname = nameof(T)

    # 1. Direct lookup in reverse table (fastest path — atomic models)
    if haskey(_XSPEC_REVERSE, tname)
        return _XSPEC_REVERSE[tname]
    end

    # 2. Must be a CompositeModel — extract left, right, and operator
    if !is_composite(model)
        throw(ArgumentError(
            "Cannot convert $(T) to XSPEC string: type not found in XSPEC_MODEL_NAMES.\n" *
            "  Only models defined in XSPECModels.jl are supported."
        ))
    end

    # CompositeModel{T,K,Operator,M1,M2}
    # type parameter 3 (index 3) is the Operator
    type_params  = T.parameters
    OperatorType = type_params[3]   # AdditionOperator, MultiplicationOperator, etc.

    is_add = OperatorType <: AdditionOperator
    # Both Multiplication and Convolution use "*" syntax in XSPEC strings
    is_mul = OperatorType <: MultiplicationOperator || OperatorType <: ConvolutionOperator

    left_model  = getfield(model, :left)
    right_model = getfield(model, :right)

    left_str  = _to_xspec(left_model,  false)
    # The right operand needs parentheses when:
    #   - The current operation is * (or conv)
    #   - AND the right sub-expression is a + (otherwise precedence is wrong)
    right_needs_parens = is_mul && is_composite(right_model) &&
                         (typeof(right_model).parameters[3] <: AdditionOperator)
    right_str = _to_xspec(right_model, right_needs_parens)

    result = if is_mul
        "$left_str*$right_str"
    else   # is_add (or any unknown op → fallback to +)
        "$left_str+$right_str"
    end

    return needs_parens ? "($result)" : result
end
