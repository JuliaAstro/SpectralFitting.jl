# Using spectral models

```@setup using_models
using SpectralFitting
using Plots
```

In this page you'll find how to use the spectral model library and how to define your own models. Using the model library is as easy as invoking or composing models from the [Model Index](@ref). For example:

```@example using_models
model = PowerLaw()
```

In the output of the REPL we see the model name, and it's two parameters, with information about those parameters, such as the current value, the associated error (10% by defaul), the minimum and maximum values, and whether the parameter is frozen or not.

!!! note

    See [`FitParam`](@ref) for full details about fitable parameters.

The parameters can be tweaked by accessing the fields

```@example using_models
model.K.value = 2.0
model.K.frozen = true
model
```

We can invoke the model on a domain in the following way

```@example using_models
domain = collect(range(0.1, 10.0, 100))
invokemodel(domain, model)
```

!!! note
    By default, models are implemented to accept a single input vector with all of the low and high bin edges, and return a flux array with the flux in each energy bin. As such, it is here the case that:
    ```julia
    length(flux) == length(energy) - 1
    ```
    Models need not be defined as such, however. See [`AbstractDataLayout`](@ref) for more.

Models can be composed together following the [Model algebra](@ref). That means to expressive a photoelectric absorption component acting on the power law we can write

```@example using_models
model2 = PhotoelectricAbsorption() * model
```

The parameters of this [`CompositeModel`](@ref) are are copied from the expression. This means we can modify the `K_1` parameter in `model2` without having to worry that we are changing `model.K`:

```@example using_models
model2.a1.K.frozen = false
model2
```

```@example using_models
model
```

Composite models have the same methods as single models. This means we can invoke a model in the same way

```@example using_models
invokemodel(domain, model2)
```

## Defining new models

To define your own model, you need to tell the package what the model parameters are and how to invoke the model. This is all done by creating a `struct` which subtypes [`AbstractSpectralModel`](@ref).

Let's create a new [`Additive`](@ref) spectral model:

```@example using_models
Base.@kwdef struct MyModel{T} <: AbstractSpectralModel{T,Additive}
    K::T = FitParam(2.0)
    p::T = FitParam(3.0)
end

# implementing a dummy add operation this function can do anything it likes, but
# must write the output into `output` and ideally should be thread safe
function SpectralFitting.invoke!(output, input, model::MyModel)
    SpectralFitting.finite_diff_kernel!(output, input) do E
        E + model.p
    end
end
```

Here we used the utility method [`SpectralFitting.finite_diff_kernel!`](@ref) to ensure the additive model is appropriately scaled across the bin width.

Note that [`Additive`](@ref) models do not need to use the normalization parameter `K` themselves. This is because when we use [`invokemodel`](@ref) these sorts of translations are automatically applied, for compatability with external models.

Our model is now ready to use
```@example using_models
model = MyModel()
```

```@example using_models
domain = collect(range(0.1, 10.0, 100))
invokemodel(domain, model)
```

!!! note
    To add new XSPEC or foreign function models, see [Wrapping new XSPEC models](@ref).

## Model abstraction

All spectral models are a sub-type of [`AbstractSpectralModel`](@ref).

```@docs
AbstractSpectralModel
SpectralFitting.invoke!
implementation
```

## Model methods

```@docs
invokemodel
invokemodel!
```

## Model algebra

Models exist as three different kinds, defined by an [`AbstractSpectralModelKind`](@ref) trait.
```@docs
AbstractSpectralModelKind
Additive
Multiplicative
Convolutional
```
## Model data availability

Many of the XSPEC implemented models use tabular data, such as FITS, and return results interpolated from these pre-calculated tables. In some cases, these table models have data files that are multiple gigabytes in size, and would be very unwieldy to ship indiscriminantly. SpectralFitting attempts to circumnavigate this bloat by downloading the model data on an _ut opus_ basis.

```@docs
SpectralFitting.download_model_data
SpectralFitting.download_all_model_data
```

Special care must be taken if new XSPEC models are wrapped to ensure the data is available. For more on this, see [Wrapping new XSPEC models](@ref).

Model data may also alternatively be copied in _by-hand_ from a HEASoft XSPEC source directory. In this case, the location to copy the data to may be determined via `joinpath(SpectralFitting.LibXSPEC_jll.artifact_dir, "spectral", "modelData")`.

## Writing XSPEC table models

SpectralFitting.jl provides functionality to create XSPEC-compatible table models and write them to FITS files following the OGIP 92-009 specification. This allows you to generate custom table models that can be used both in SpectralFitting.jl and in XSPEC.

```@docs
write_table_model
```

### Example: Creating a custom table model

Here's an example of creating a broken power law table model:

```julia
using SpectralFitting

# Define energy grid (500 bins, logarithmically spaced)
energy_bins = exp.(range(log(0.1), log(20.0), length=501))

# Define parameter grids
index1_grid = collect(range(1.5, 2.5, step=0.1))  # 11 values (lower energy index)
break_energy_grid = exp.(range(log(0.5), log(10.0), length=20))  # 20 values, log-spaced
index2_grid = collect(range(1.5, 2.5, step=0.1))  # 11 values (upper energy index)
param_grids = (index1_grid, break_energy_grid, index2_grid)

# Calculate spectra for each parameter combination
n_energy_bins = length(energy_bins) - 1
n_spectra = length(index1_grid) * length(break_energy_grid) * length(index2_grid)
spectra = zeros(Float64, n_energy_bins, n_spectra)

# Compute mid-point energies and bin widths
E_mid = (energy_bins[1:end-1] .+ energy_bins[2:end]) ./ 2
dE = energy_bins[2:end] .- energy_bins[1:end-1]

# Fill spectra array (first parameter varies fastest)
# Loop structure: index1 (fastest) -> E_break (middle) -> index2 (slowest)
# This matches XSPEC parameter order: PhoIndex1, BreakE, PhoIndex2
idx = 1
for index2 in index2_grid
    for E_break in break_energy_grid
        for index1 in index1_grid
            # Broken power law: E^(-index1) at low energies, E^(-index2) at high energies
            # with continuity at the break energy
            norm_factor = E_break^(index2 - index1)
            # Evaluate at bin center and integrate over bin width
            spectra[:, idx] = [(E < E_break ? E^(-index1) : norm_factor * E^(-index2)) * width
                               for (E, width) in zip(E_mid, dE)]
            idx += 1
        end
    end
end

# Specify interpolation method for each parameter
# 0 = linear interpolation, 1 = logarithmic interpolation
param_method = [Int32(0), Int32(1), Int32(0)]  # log interpolation for break energy

# Specify interpolation method for each parameter
# 0 = linear interpolation, 1 = logarithmic interpolation
param_method = [Int32(0), Int32(1), Int32(0)]  # log interpolation for break energy

# Write to FITS file
write_table_model(
    "bknpowerlaw_table.fits",
    energy_bins,
    param_grids,
    spectra;
    model_name = "BKNPOW",
    param_names = ["PhoIndex1", "BreakE", "PhoIndex2"],
    param_units = ["", "keV", ""],
    model_units = "photons/cm^2/s/keV",
    additive = true,
    param_method = param_method
)
```

The resulting FITS file can then be loaded as a table model:

```julia
tmd = TableModelData(Val(3), "bknpowerlaw_table.fits")
```

