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

"""
    write_table_model(
        path::String,
        energy_bins::AbstractVector{T},
        param_grids::NTuple{N,Vector{T}},
        spectra_array::Array{T};
        model_name::String = "TABLE",
        model_units::String = "photons/cm^2/s",
        redshift::Bool = false,
        additive::Bool = true,
        param_names::Vector{String} = String["par\$i" for i in 1:N],
        param_units::Vector{String} = fill("", N),
        param_method::Vector{Int32} = fill(Int32(0), N),
        lowE_units::String = "keV",
        highE_units::String = "keV"
    ) where {T,N}

Write an XSPEC-compatible table model to a FITS file following the OGIP 92-009 specification.

This is the inverse operation of reading a table model with [`TableModelData`](@ref).

# Arguments
- `path::String`: Output FITS file path
- `energy_bins::AbstractVector{T}`: Energy bin edges. Should have length `n_energy_bins`.
  The bins are assumed contiguous (i.e., `energy_bins[i+1] = E_hi[i]`).
- `param_grids::NTuple{N,Vector{T}}`: Tuple of parameter grid axes. Each vector contains
  the unique values for that parameter dimension.
- `spectra_array::Array{T}`: Array of spectra with shape `(n_energy_bins-1, n_spectra)` where
  `n_spectra = prod(length.(param_grids))`. The spectra should be ordered such that the
  first parameter varies fastest (C-order/row-major).

# Keywords
- `model_name::String = "TABLE"`: Name of the table model
- `model_units::String = "photons/cm^2/s"`: Units of the spectral values
- `redshift::Bool = false`: Whether this is a redshift table model
- `additive::Bool = true`: Whether this is an additive (true) or multiplicative (false) model
- `param_names::Vector{String}`: Names for each parameter (length N)
- `param_units::Vector{String}`: Units for each parameter (length N)
- `param_method::Vector{Int32}`: Interpolation method for each parameter (length N). 
  Use 0 for linear interpolation, 1 for logarithmic interpolation. Defaults to linear (0) for all parameters.
- `lowE_units::String = "keV"`: Units for low energy bin edges
- `highE_units::String = "keV"`: Units for high energy bin edges

# Example
```julia
# Create a simple 2-parameter table model
energy_bins = collect(range(0.1, 10.0, length=101))  # 100 energy bins
param1_grid = [1.0, 2.0, 3.0]
param2_grid = [0.5, 1.0, 1.5, 2.0]
param_grids = (param1_grid, param2_grid)

# Generate spectra for each parameter combination
n_energy_bins = length(energy_bins) - 1
n_spectra = length(param1_grid) * length(param2_grid)
spectra = zeros(n_energy_bins, n_spectra)

# Fill with model values (example: simple power law-like model)
idx = 1
for p2 in param2_grid
    for p1 in param1_grid
        E_mid = (energy_bins[1:end-1] .+ energy_bins[2:end]) ./ 2
        spectra[:, idx] = p1 .* E_mid.^(-p2)
        idx += 1
    end
end

# Write with default linear interpolation for all parameters
write_table_model(
    "my_table_model.fits",
    energy_bins,
    param_grids,
    spectra;
    model_name = "MYTABLE",
    param_names = ["index", "xi"]
)

# Or specify interpolation method: linear for index, logarithmic for xi
write_table_model(
    "my_table_model.fits",
    energy_bins,
    param_grids,
    spectra;
    model_name = "MYTABLE",
    param_names = ["index", "xi"],
    param_method = Int32[0, 1]  # 0=linear, 1=logarithmic
)
```
"""
function write_table_model(
    path::String,
    energy_bins::AbstractVector{T},
    param_grids::NTuple{N,Vector{T}},
    spectra_array::Array{T};
    model_name::String = "TABLE",
    model_units::String = "photons/cm^2/s",
    redshift::Bool = false,
    additive::Bool = true,
    param_names::Vector{String} = String["par$i" for i in 1:N],
    param_units::Vector{String} = fill("", N),
    param_method::Vector{Int32} = fill(Int32(0), N),
    lowE_units::String = "keV",
    highE_units::String = "keV"
) where {T,N}
    # Validate inputs
    n_energy_bins = length(energy_bins) - 1
    @assert n_energy_bins > 0 "Energy bins must have at least 2 edges"
    @assert issorted(energy_bins) "Energy bins must be sorted"
    
    n_spectra = prod(length.(param_grids))
    @assert size(spectra_array, 1) == n_energy_bins "Spectra array first dimension must match number of energy bins"
    @assert size(spectra_array, 2) == n_spectra "Spectra array second dimension must match total number of parameter combinations"
    @assert length(param_names) == N "Number of parameter names must match number of parameter grids"
    @assert length(param_units) == N "Number of parameter units must match number of parameter grids"
    @assert length(param_method) == N "Number of parameter methods must match number of parameter grids"
    @assert all(m -> m == 0 || m == 1, param_method) "Parameter method values must be 0 (linear) or 1 (logarithmic)"
    
    # Create FITS file
    FITS(path, "w") do f   
        # HDU 1: Primary (empty but with key headers)
        write(f, Int[])
        
        # Add key headers to primary HDU
        primary = f[1]
        write_key(primary, "HDUCLASS", "OGIP", "Conforms to OGIP standard")
        write_key(primary, "HDUCLAS1", "XSPEC TABLE MODEL", "Model spectra for XSPEC")
        write_key(primary, "MODLNAME", model_name, "Model name")
        write_key(primary, "MODLUNIT", model_units, "Model units")
        write_key(primary, "REDSHIFT", redshift, "If true then redshift will be a parameter")
        write_key(primary, "ADDMODEL", additive, "If true then this is additive table model")
        write_key(primary, "HDUVERS", "1.0.0", "Version of format")
        write_key(primary, "NINTPARM", Int32(N), "Number of interpolation parameters")
        write_key(primary, "NADDPARM", Int32(0), "Number of additional parameters")
        
        # HDU 2: PARAMETERS extension
        param_names_arr = collect(param_names)
        param_initial = [first(pg) for pg in param_grids]
        param_bottom = [minimum(pg) for pg in param_grids]
        param_top = [maximum(pg) for pg in param_grids]
        param_numbvals = [length(pg) for pg in param_grids]
        param_delta = [length(pg) > 1 ? (pg[2] - pg[1]) : zero(T) for pg in param_grids]
        
        # Build units dictionary - one per column
        param_col_units = Dict{String,String}()
        for (i, name) in enumerate(param_names)
            param_col_units["INITIAL"] = get(param_units, 1, "")
            param_col_units["DELTA"] = get(param_units, 1, "")
            param_col_units["MINIMUM"] = get(param_units, 1, "")
            param_col_units["BOTTOM"] = get(param_units, 1, "")
            param_col_units["TOP"] = get(param_units, 1, "")
            param_col_units["MAXIMUM"] = get(param_units, 1, "")
            param_col_units["VALUE"] = get(param_units, 1, "")
        end
        
        # Convert param_grids to Vector{Vector{T}} for VALUE column (variable-length arrays)
        param_values = [collect(pg) for pg in param_grids]
        
        write(f, ["NAME", "METHOD", "INITIAL", "DELTA", "MINIMUM", "BOTTOM", "TOP", "MAXIMUM", "NUMBVALS", "VALUE"],
              [param_names_arr, param_method, param_initial, param_delta, 
               param_bottom, param_bottom, param_top, param_top, 
               Int32.(param_numbvals), param_values];
              units = param_col_units,
              varcols = ["VALUE"],
              name = "PARAMETERS")
        
        # Add headers after table creation
        hdu = f[end]
        write_key(hdu, "HDUCLASS", "OGIP", "Conforms to OGIP standard")
        write_key(hdu, "HDUCLAS1", "XSPEC TABLE MODEL", "Model spectra for XSPEC")
        write_key(hdu, "HDUCLAS2", "PARAMETERS", "Extension containing parameter info")
        write_key(hdu, "HDUVERS", "1.0.0", "Version of format")
        write_key(hdu, "NINTPARM", Int32(N), "Number of interpolation parameters")
        write_key(hdu, "NADDPARM", Int32(0), "Number of additional parameters")
        write_key(hdu, "REDSHIFT", redshift, "If true then redshift will be a par")
        write_key(hdu, "ADDMODEL", additive, "If true then this is additive table model")
        
        # HDU 3: ENERGIES extension  
        energy_lo = energy_bins[1:end-1]
        energy_hi = energy_bins[2:end]
        
        write(f, ["ENERG_LO", "ENERG_HI"],
              [collect(energy_lo), collect(energy_hi)];
              units = Dict("ENERG_LO" => lowE_units, "ENERG_HI" => highE_units),
              name = "ENERGIES")
        
        # Add headers after table creation
        hdu = f[end]
        write_key(hdu, "HDUCLASS", "OGIP", "Conforms to OGIP standard")
        write_key(hdu, "HDUCLAS1", "XSPEC TABLE MODEL", "Model spectra for XSPEC")
        write_key(hdu, "HDUCLAS2", "ENERGIES", "Extension containing energies")
        write_key(hdu, "HDUVERS", "1.0.0", "Version of format")
        
        # HDU 4: SPECTRA extension
        # Create the parameter value mesh and reorder spectra for XSPEC
        # XSPEC expects last parameter to vary fastest, but input has first parameter fastest
        # So we need to reverse the parameter order
        param_mesh = zeros(T, N, n_spectra)
        reordered_spectra = zeros(T, size(spectra_array))
        
        # Create mapping from input order (first param fastest) to XSPEC order (last param fastest)
        input_indices = CartesianIndices(Tuple(length.(param_grids)))
        reversed_lengths = Tuple(reverse([length(pg) for pg in param_grids]))
        output_indices = CartesianIndices(reversed_lengths)
        
        for (input_idx, input_cart) in enumerate(input_indices)
            # Reverse the multi-dimensional index for output
            reversed_cart = reverse(Tuple(input_cart))
            output_idx = LinearIndices(reversed_lengths)[CartesianIndex(reversed_cart)]
            
            # Copy spectrum to new position
            reordered_spectra[:, output_idx] = spectra_array[:, input_idx]
            
            # Store parameters in reversed order (last parameter in first row, etc.)
            for (param_idx, grid_idx) in enumerate(Tuple(input_cart))
                param_mesh[N - param_idx + 1, output_idx] = param_grids[param_idx][grid_idx]
            end
        end
        
        write(f, ["PARAMVAL", "INTPSPEC"],
              [param_mesh, reordered_spectra];
              units = Dict("INTPSPEC" => model_units),
              name = "SPECTRA")
        
        # Add headers after table creation
        hdu = f[end]
        write_key(hdu, "HDUCLASS", "OGIP", "Conforms to OGIP standard")
        write_key(hdu, "HDUCLAS1", "XSPEC TABLE MODEL", "Model spectra for XSPEC")
        write_key(hdu, "HDUCLAS2", "MODEL SPECTRA", "Extension containing model spectra")
        write_key(hdu, "HDUVERS", "1.0.0", "Version of format")
        write_key(hdu, "MODLNAME", model_name, "Model name")
        write_key(hdu, "MODLUNIT", model_units, "Model units")
        write_key(hdu, "REDSHIFT", redshift, "If true then redshift will be a parameter")
        write_key(hdu, "ADDMODEL", additive, "If true then this is additive table model")
    end
    
    return path
end

export AbstractTableModel, TableModelData, TableGridData, TableModelInterpolation,
       write_table_model
