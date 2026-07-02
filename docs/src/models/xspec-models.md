```@meta
CurrentModule = XSPECModels
```

# XSPEC models

SpectralFitting can vendor the XSPEC model library through the utility package XSPECModels.jl. You can install this package directly from the University of Bristol's [AstroRegistry](https://github.com/astro-group-bristol/AstroRegistry):

```julia
pkg> registry add https://github.com/astro-group-bristol/AstroRegistry
pkg> add XSPECModels
```

Using the package is straight forward once installed:
```julia
using SpectralFitting, XSPECModels
```

!!! note
    The convention is that models that have are imported from XSPEC or have XSPEC ABI are prefixed with `XS_` in their name. For example, the XSPEC equivalent of [`PowerLaw`](@ref) is [`XS_PowerLaw`](@ref).


The XSPEC models are wrapped into SpectralFitting models using the [`XSPECModels.@xspecmodel`] macro. If a model you require is not already wrapped, this macro will make that easy to do. Please consider upstreaming your wrapper via a PR to the SpectralFitting GitHub repository.

```@docs
@xspecmodel
```

XSPEC models frequently have tabular data dependencies, without which the models fail to invoke (see [Model data availability](@ref)). If the data files are known but not present, the XSPEC models will throw an error with instructions for downloading the missing data. If the data files are unknown, Julia may crash catastrophically. If this is the case, often a single line will be printed with the LibXSPEC error, specifying the name of the missing source file. This can be registered as a data dependency of a model using [`SpectralFitting.register_model_data`](@ref).

The first time any XSPEC model is invoked, SpectralFitting checks to see whether requisite data is needed, and whether the data is downloaded. Subsequent calls will hit a lookup cache instead to avoid run-time costs of performing this check.

## Additive


```@docs
XS_PowerLaw
XS_BlackBody
XS_BremsStrahlung
XS_Gaussian
XS_Laor
XS_DiskLine
XS_PhotoelectricAbsorption
XS_WarmAbsorption
XS_CalculateFlux
XS_KerrDisk
XS_KyrLine
```
