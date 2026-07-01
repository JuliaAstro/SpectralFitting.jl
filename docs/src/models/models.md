# Model index

Models wrapped from XSPEC implementations are prefixed with `XS_*`, whereas pure-Julia models are simply named, e.g. [`XS_PowerLaw`](@ref) in XSPEC vs [`PowerLaw`](@ref) in Julia.

The available models are
```@index
Pages = ["models.md"]
Order = [:type]
```

## Julia models

### Additive

```@docs
BlackBody
BremsStrahlung
DeltaLine
GaussianLine
PowerLaw
```

### Multiplicative
```@docs
PhotoelectricAbsorption
```

### Convolutional
```@docs
Log10Flux
```

## Utility

```@docs
SpectralFitting.register_model_data
SpectralFitting.finite_diff_kernel!
wrap_model_as_objective
```

## Generating model fingerprints

To generate the unicode plot to add as a fingerprint, we use a simple function:

```@example fingerprints
using SpectralFitting, UnicodePlots

function plotmodel(energy, model)
    flux = invokemodel(energy, model)
    lineplot(
        energy[1:end-1],
        flux,
        title=String(Base.typename(typeof(model)).name),
        xlabel="E (keV)",
        canvas=DotCanvas
    )
end

# e.g. for XS_PowerLaw()
energy = collect(range(0.1, 20.0, 100))
plotmodel(energy, PowerLaw())
```
