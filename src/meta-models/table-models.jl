"""
    abstract type AbstractTableModel{T,K} <: AbstractSpectralModel{T,K} end

Abstract type representing table models, i.e. those models that interpolate or
load data from a table.

First field in the struct **must be** `table`. See
[`PhotoelectricAbsorption`](@ref) for an example implementation.
"""
abstract type AbstractTableModel{T,K} <: AbstractSpectralModel{T,K} end

"""
    Base.copy(m::AbstractTableModel)

Create a copy of an [`AbstractTableModel`](@ref). This will copy all fields except
the `table` field, which is assumed to be a constant set of values that can be
shared by multiple copies.

When this is not the case, the user should redefine `Base.copy` for their particular
table model to copy the table as needed.
"""
function Base.copy(m::AbstractTableModel)
    typeof(m)(m.table, (copy(getproperty(m, f)) for f in fieldnames(typeof(m))[2:end])...)
end

abstract type AbstractCachedModel{T,K} <: AbstractSpectralModel{T,K} end

# some utilities for interacting with XSPEC-compatible table models

"""
    TableGridData

Used to contain table grid data in the interpolation scheme.
"""
struct TableGridData{V}
    "Interpolated values or values to be interpolated, depending on calling semantics."
    values::V
end

MultiLinearInterpolations.restructure(::TableGridData, vs::AbstractVector) =
    TableGridData(vs)

"""
    TableModelData

A container for data read in from an XSPEC-style `Xtable` table model. This
structure currently assumes a number of things:
- The energy bins are strictly contiguous (`E_lo[i+1] = E_hi[i]`)
- The parameter grid is rectilinear

$(FIELDS)
"""
struct TableModelData{T,N}
    "Energy bins used for all of the grid entries."
    energy_bins::Vector{T}
    "The parameter axes, with the range of tabulated values."
    params::NTuple{N,Vector{T}}
    "All grids, laid out in the same way as the parameter axes."
    grids::Array{TableGridData{Vector{T}},N}
end

function Base.show(
    io::IO,
    ::MIME"text/plain",
    @nospecialize(data::TableModelData{T,N})
) where {T,N}
    print(
        io,
        "TableModelData{$T}[N=$(N),E_bins=$(length(data.energy_bins)),grid_size=$(length(data.grids))]",
    )
end

"""
    TableModelInterpolation
    TableModelInterpolation(tmd::TableModelData)

Wraps a [`TableModelData`](@ref) and augments it with a multi-linear
interpolation cache.

This can then be used with `interpolate_table!` to interpolate the table data.
"""
struct TableModelInterpolation{T,N}
    data::TableModelData{T,N}
    interpolator::MultilinearInterpolator{N,T}
end

function Base.show(
    io::IO,
    ::MIME"text/plain",
    interp::TableModelInterpolation{T,N},
) where {T,N}
    print(io, "TableModelInterpolation($T,$N)")
end

TableModelInterpolation(tmd::TableModelData{T,N}) where {T,N} =
    TableModelInterpolation(tmd, MultilinearInterpolator{N}(tmd.grids; T = T))

function interpolate_table!(tmi::TableModelInterpolation{T}, parameters::Vararg) where {T}
    v = MultiLinearInterpolations.interpolate!(
        tmi.interpolator,
        tmi.data.params,
        tmi.data.grids,
        reverse(convert.(T, parameters)),
    )
    v.values
end

"""
    TableModelData(::Val{N}, path::String; T::Type = Float64)

Unpack an XSPEC table model into a [`TableModelData`](@ref) for use in wrapping
table models. The first argument must be compile time known, and be the number
of parameter the model has. SpectralFitting.jl needs to know this ahead-of-time
as part of it's fitting optimisations. The `path` is the path to the
corresponding fits file, and the `T` keyword argument is what value type to cast
the data into.

Many models will use `Float32`, which the caller may want to consider, if the
tables are very large.
"""
function TableModelData(::Val{N}, path::String; T::Type = Float64) where {N}
    f = fits(path)

    # this assumes energy is always contiguously binned, probably valid enough
    # for XSPEC table models
    energy_bins::Vector{T} = convert.(T, read(f[3], "ENERG_LO"))
    @assert issorted(energy_bins) "Energy bins are not in a linear order"
    _energy_high = convert.(T, read(f[3], "ENERG_HI"))
    @views @assert all(isapprox.(_energy_high[1:(end-1)], energy_bins[2:end])) "Energy bins are not contiguously binned (i.e. `e_low[i + 1] != e_high[i]`)."
    push!(energy_bins, last(_energy_high))

    # currently assuming the actual underlying parameter grid is
    # rectilinear. if it isn't this parsing will fail.
    parameter_meshgrid::Matrix{T} = convert.(T, read(f[4], "PARAMVAL"))
    parameter_axes = map(1:size(parameter_meshgrid, 1)) do i
        unique!(parameter_meshgrid[i, :])
    end
    reverse!(parameter_axes)

    _data::Matrix{T} = convert.(T, read(f[4], "INTPSPEC"))
    _expected_size = reduce(*, length(i) for i in parameter_axes)
    # this assert should hopefully be good enough to diagnose errors related to
    # things not being rectilinear
    @assert size(_data, 2) == _expected_size "Parameter axes are unlikely recti-linear, else bad number of spectra in `INTPSPEC`: size(_data, 2) = $(size(_data, 2)) != $(_expected_size)"
    _grid = map(1:size(_data, 2)) do i
        TableGridData(_data[:, i])
    end

    param_tuple = ((parameter_axes[i] for i = 1:N)...,)
    grid = reshape(_grid, length.(param_tuple))

    TableModelData{T,N}(energy_bins, param_tuple, grid)
end

export AbstractTableModel, TableModelData, TableGridData, TableModelInterpolation
