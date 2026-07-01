# Parameters

One of the core interface that SpectralFitting.jl provides are the model
abstractions, and a method for controlling their _parameters_. Models are
defined with a generic field type, so that they can be converted to primitive
types (e.g. Float64 or Float32) during evaluation. But whilst defining a model
or a problem, the field type of all parameters is [`FitParam`](@ref).

```@docs 
FitParam
```

## Parameter binding

When performing a fit, it is desirable to **bind** certain parameters together. This ensures that they will have the same value; for example, if you were fitting two simultaneous datasets with two [`PowerLaw`](@ref) models, you may want to have different normalisations of the model components, but enforce the power law index to be the same. To achieve this, SpectralFitting has the [`bind!`](@ref) function that applies to your [`FittingProblem`](@ref).

```@docs
bind!
bindall!
```

!!! note
    Bindings are treated not as specific to the model but specific to the [`FittingProblem`](@ref). This is because you may want to use the same model for multiple different datasets, and have slightly different binding requirements for each one (e.g. depending on the instruments you are using). If you do need the same binding applied to two different problems, you can do that with
    ```julia
    append!(prob1.bindings, prob2.bindings)
    ```
    Caution however, this will only make sense if you are using precisely the same model in both problems.


Let's try it out. We'll generate some arbitrary powerlaw spectra with different normalisations and fit them simultaneously.
```@example bind
using SpectralFitting, Plots

energy = collect(range(0.1, 10.0, 100))

# two different models with different normalisations
model1 = PowerLaw(K = FitParam(100.0), a = FitParam(1.2))
model2 = PowerLaw(K = FitParam(300.0), a = FitParam(2.0))

data1 = simulate(energy, model1, var = 1e-3, seed = 42)
data2 = simulate(energy, model2, var = 1e-3, seed = 42)

plot(data1, xscale = :log10, yscale = :log10)
plot!(data2, xscale = :log10, yscale = :log10)
```

Now we want to fit a single powerlaw model to both of these spectra simultaneously, but with the powerlaw index fixed to be the same in both models.
```@example bind
model = PowerLaw()
prob = FittingProblem(model => data1, model => data2)

bindall!(prob, :a)
prob
```

We can get a better look at our model configuration by using the [`details`](@ref) method:
```@example bind
details(prob)
```

In this printout we see that the `a` parameter of `Model 2` is bound to the `a` parameter of `Model 1`.

```@example bind
result = fit(prob, LevenbergMarquadt())

plot(data1, xscale = :log10, yscale = :log10)
plot!(data2, xscale = :log10, yscale = :log10)
plot!(result[1])
plot!(result[2])
```

Note that this fit is bad, because the underlying data have different power law indices, but our fit is required to enforce the models to have the same value. If we release this requirement, the fit will be much better, but the models will be entirely independent.

```@example bind
prob = FittingProblem(model => data1, model => data2)

result = SpectralFitting.fit(prob, LevenbergMarquadt())

plot(data1, xscale = :log10, yscale = :log10)
plot!(data2, xscale = :log10, yscale = :log10)
plot!(result[1])
plot!(result[2])
```

## Parameter patching

!!! warning

    The parameter patching interface is still experimental and will likely change in subsesquent versions of SpectralFitting.jl

Sometimes simply linking values is not sufficient, and you need to express a complex relationship between parameters. This is where [`ParameterPatch`](@ref), a type of [`AbstractModelWrapper`](@ref) is useful.

Consider the following

```@example bind
using SpectralFitting, Plots

real_model = GaussianLine(K = FitParam(2.0)) + GaussianLine(μ = FitParam(2.0), K = FitParam(4.0))
energy = collect(range(0.1, 10.0, 100))

data = simulate(energy, real_model; var = 4e-4, seed = 42)
plot(data)
```

We can fit this simple dataset quite easily, but suppose we wanted to constrain the solution to have the normalisation of one model component be exactly twice that of the other. To achieve this, we can use a patch:

```@example bind
function my_patch!(p)
    p.a2.K = p.a1.K * 2
end

# reset the parameters so the fit starts "fresh"
model = GaussianLine() + GaussianLine(μ = FitParam(3.0))

# Wrap the model with a parameter patch
patched_model = ParameterPatch(model; patch = my_patch!)
# be sure to freeze any parameter you are planning to overwrite in a patch
patched_model.a2.K.frozen = true

patched_model
```

This can then be fit as usual:

```@example bind
prob = FittingProblem(patched_model => data)
result = fit(prob, LevenbergMarquadt())

plot!(result)
```

To apply the result with a parameter patch back on a model, use [`update_model!`](@ref) or [`apply_patch!`](@ref)

```@example bind
update_model!(patched_model, result)
patched_model
```