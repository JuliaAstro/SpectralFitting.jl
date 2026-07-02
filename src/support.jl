@enumx ErrorStatistics begin
    Numeric
    Poisson
    Gaussian
    Unknown
end

"""
    abstract type AbstractDataLayout end

The data layout primarily concerns the relationship between the objective and
the domain. It is used to work out whether a model and a dataset are fittable,
and if not, whether a translation in the output of the model to the domain of
the model is possible.

The following methods may be used to interrogate support:
- [`preferred_support`](@ref) for inferring the preferred support of a model when multiple supports are possible.
- [`common_support`](@ref) to obtain the common support of two structures

The following method is also used to define the support of a model or dataset:
- [`supports`](@ref)

For cases where unit information needs to be propagated, an `AbstractDataLayout` can also be used to ensure the units are compatible. To query the units of a layout, use
- [`support_units`](@ref)
"""
abstract type AbstractDataLayout end

"""
    support_units(x)

Return the units of a particular layout. If this method returns `nothing`,
assume the layout does not care about the units and handle that information
appropriately (throw an error or set defaults).
"""
support_units(::AbstractDataLayout) = nothing

"""
    with_units(::AbstractDataLayout, units)

Remake the [`AbstractDataLayout`](@ref) with the desired units. This may be a no-op
if the layout does not care about units, see [`support_units`](@ref).
"""
with_units(layout::AbstractDataLayout, units) = layout

default_layout(data) =
    with_units(preferred_support(data), preferred_units(data, ChiSquared()))
default_layout(model, data) =
    with_units(preferred_support(model), preferred_units(data, ChiSquared()))

"""
    struct OneToOne <: AbstractDataLayout end

Indicates there is a one-to-one (injective) correspondence between each input
value and each output value. That is to say
```julia
length(objective) == length(domain)
```
"""
struct OneToOne{U} <: AbstractDataLayout
    units::U
end
OneToOne() = OneToOne(nothing)
support_units(l::OneToOne) = l.units
with_units(::OneToOne, units) = OneToOne(units)

"""
    struct ContiguouslyBinned <: AbstractDataLayout end

Contiguously binned data layout means that the domain describes high and low
bins, with the objective being the value in that bin. This means
```julia
length(objective) + 1== length(domain)
```
Note that the _contiguous_ qualifier is to mean there is no gaps in the bins, and that
```math
ΔE_i = E_{i+1} - E_{i}
```

"""
struct ContiguouslyBinned{U} <: AbstractDataLayout
    units::U
end
ContiguouslyBinned() = ContiguouslyBinned(nothing)
support_units(l::ContiguouslyBinned) = l.units
with_units(::ContiguouslyBinned, units) = ContiguouslyBinned(units)

const DEFAULT_SUPPORT_ORDERING = (ContiguouslyBinned(), OneToOne())

"""
    preferred_support(x)

Get the preferred [`AbstractDataLayout`](@ref) of `x`. If multiple supports are available,
the `DEFAULT_SUPPORT_ORDERING` is followed:

```
DEFAULT_SUPPORT_ORDERING = $(DEFAULT_SUPPORT_ORDERING)
```
"""
function preferred_support(x)
    for layout in DEFAULT_SUPPORT_ORDERING
        if supports(layout, x)
            return layout
        end
    end
    error("No preferred support for $(typeof(x))")
end

"""
    common_support(x, y)

Find the common [`AbstractDataLayout`](@ref) of `x` and `y`, following the ordering of
[`preferred_support`](@ref).
"""
function common_support(x, y)
    # todo: check if can be compile time eval'd else expand for loop or @generated

    # order of preference is important
    # can normally trivially fallback from one-to-one to contiguous bins to regular bins
    for layout in DEFAULT_SUPPORT_ORDERING
        if supports(layout, x) && supports(layout, y)
            return layout
        end
    end
    error("No common support between $(typeof(x)) and $(typeof(y)).")
end

function _support_reducer(x::OneToOne, y)
    if supports(x, y)
        return x
    else
        error("No common support!!")
    end
end
function _support_reducer(x::ContiguouslyBinned, y)
    if supports(x, y)
        return x
    else
        _support_reducer(OneToOne(), y)
    end
end
_support_reducer(x, y) = common_support(x, y)

common_support(args::Vararg) = reduce(_support_reducer, args)

"""
    supports(x::Type)

Used to define whether a given type has support for a specific
[`AbstractDataLayout`](@ref). Should return a tuple of the supported layouts. This
method should be implemented to express new support, not the query method.

To query, there is

    supports(layout::AbstractDataLayout, x)::Bool

# Example

```julia
supports(::Type{typeof(x)}) = (OneToOne(),)
@assert supports(ContiguouslyBinned(), x) == false
```
"""
supports(layout::AbstractDataLayout, x)::Bool = layout in supports(x)
supports(::T) where {T} = supports(T)
supports(::Type) = ()

export AbstractDataLayout,
    common_support, ContiguouslyBinned, OneToOne, preferred_support, supports, support_units
