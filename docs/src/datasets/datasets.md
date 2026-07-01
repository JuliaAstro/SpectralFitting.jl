# Datasets

SpectralFitting.jl supports a wide variety of datasets, and makes it easy to wrap your own.

For spectral fitting specifics, the main dataset type is

```@docs
SpectralData
```

## Dataset abstraction

Datasets must define a small API to make fitting possible. The picture to have in mind when considering the different domains is as follows: the model is trying to predict the objective. It does so by taking in input domain and maps it to some output domain.

That means [`make_output_domain`](@ref) and [`make_objective_domain`](@ref) correspond to the $(X,Y)$ values of the data that the model is trying to fit, whilst the model is evaluated on the [`make_model_domain`](@ref), which need not be the same as the output domain.

In other cases, the [`objective_transformer`](@ref) acts to transform the output of the model onto the output domain. 

Mathematically, expressing the output domain $X$, the model domain $D$, the model output $M(D)$ and objective $S$, along with the transformer as $T$, then the relationship between the different domains is

```math
\hat{S} = T \times M(D),
```

Both $\hat{S}$ and $S$ are defined over $X$. The various fitting operations try to find model parameters that make $\hat{S}$ and $S$ as close as possible.

```@docs
AbstractDataset
make_objective_variance
make_objective
make_domain_variance
make_model_domain
make_output_domain
```

## Underlying data layouts

```@docs
AbstractDataLayout
OneToOne
ContiguouslyBinned
common_support
preferred_support
supports
```