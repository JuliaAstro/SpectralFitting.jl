# Transitioning from XSPEC

This page is for astronomers familiar with
[XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/) who want to use
SpectralFitting.jl. The two packages share the same underlying model library
(via [LibXSPEC_jll](https://github.com/astro-group-bristol/LibXSPEC_jll.jl)),
but SpectralFitting.jl expresses models as Julia types rather than strings.

## Model string translation

If you have an existing XSPEC model string, you can translate it directly into
a SpectralFitting.jl model using the `xspec"..."` string macro provided by
XSPECModels.jl:

```julia
using SpectralFitting, XSPECModels

# Classic absorbed power-law
model = xspec"phabs*powerlaw"

# Thermal plasma with neutral hydrogen absorption  
model = xspec"tbabs*powerlaw"

# Multi-component: absorption acting on sum of additive models
model = xspec"phabs*(powerlaw+bbody)"

# Three-component soft-state model
model = xspec"tbabs*(diskbb+powerlaw+gaussian)"
```

The string macro is **case-insensitive** and **whitespace-tolerant**:

```julia
xspec"PHABS * PowerLaw"   # same as xspec"phabs*powerlaw"
```

!!! note "Operator precedence"
    `*` binds more tightly than `+`, matching XSPEC conventions:
    ```julia
    # These are different models!
    xspec"phabs*powerlaw+bbody"    # → (phabs*powerlaw) + bbody  — bbody is unabsorbed
    xspec"phabs*(powerlaw+bbody)"  # → phabs * (powerlaw+bbody)  — both are absorbed
    ```
    This is astrophysically significant: only use parentheses when absorption
    should apply to all additive components.

For runtime strings (e.g. read from a file), use `parse_xspec_model_string`:

```julia
model_str = "phabs*powerlaw"   # could come from a config file
expr  = parse_xspec_model_string(model_str)
model = eval(expr)
```

## Translating a model back to an XSPEC string

Given any SpectralFitting.jl model built from XSPECModels types, you can
recover the XSPEC string representation with `xspec_model_string`:

```julia
model = XS_PhotoelectricAbsorption() * (XS_PowerLaw() + XS_BlackBody())
xspec_model_string(model)   # → "phabs*(powerlaw+bbody)"
```

This is useful for logging, reproducibility, or interoperability with
scripts that call XSPEC directly.

## XSPEC name reference

The full mapping of XSPEC model names to Julia types is in `XSPEC_MODEL_NAMES`:

```julia
XSPEC_MODEL_NAMES["tbabs"]    # :XS_NeutralHydrogenAbsorption
XSPEC_MODEL_NAMES["phabs"]    # :XS_PhotoelectricAbsorption
XSPEC_MODEL_NAMES["powerlaw"] # :XS_PowerLaw
XSPEC_MODEL_NAMES["diskbb"]   # :XS_DiskBlackBody
```

Common XSPEC abbreviations are also supported: `"po"` for `"powerlaw"`,
`"gaus"` for `"gaussian"`, `"bb"` for `"bbody"`.

## Comparison: XSPEC syntax vs SpectralFitting.jl

| XSPEC | SpectralFitting.jl |
|:------|:-------------------|
| `model phabs*powerlaw` | `xspec"phabs*powerlaw"` or `XS_PhotoelectricAbsorption() * XS_PowerLaw()` |
| `model tbabs*(apec+powerlaw)` | `xspec"tbabs*(powerlaw+bbody)"` |
| `newpar 1 0.1` | `model.ηH_1.value = 0.1` |
| `freeze 1` | `model.ηH_1.frozen = true` |
| `thaw 2` | `model.a_1.frozen = false` |

## API reference

```@docs
parse_xspec_model_string
xspec_model_string
XSPEC_MODEL_NAMES
@xspec_str
```
