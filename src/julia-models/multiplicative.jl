"""
    PhotoelectricAbsorption

$(FIELDS)
"""
struct PhotoelectricAbsorption{D,T} <: AbstractTableModel{T,Multiplicative}
    table::D
    "Equivalent hydrogen column (units of 10²² atoms per cm⁻²)."
    ηH::T
end
function PhotoelectricAbsorption(ηH::T) where {T}
    data = get_model_data(PhotoelectricAbsorption)
    E::Vector{Float64} = data[1]["E"]
    σ::Vector{Float64} = data[1]["σ"]
    table = linear_interpolation(E, σ, extrapolation_bc = Line())
    PhotoelectricAbsorption{typeof(table),T}(table, ηH)
end
PhotoelectricAbsorption(; ηH = FitParam(1.0)) = PhotoelectricAbsorption(ηH)
register_model_data(PhotoelectricAbsorption, "cross_sections_phabs_angr.jld")
@inline @fastmath function invoke!(flux, energy, model::PhotoelectricAbsorption)
    let ηH = model.ηH, table = model.table
        E = @views energy[1:(end-1)]
        @. flux = exp(-ηH * table(E))
    end
end

"""
    Constant

$(FIELDS)
"""
struct Constant{T} <: AbstractSpectralModel{T,Multiplicative}
    "Constant value."
    value::T
end
Constant(; value = FitParam(1.0; frozen = true)) = Constant(value)
@fastmath function invoke!(flux, energy, model::Constant)
    flux .= model.value
end

export PhotoelectricAbsorption, Constant
