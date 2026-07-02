"""
    Log10Flux

Used to measure the (log) flux of the models it is applied to. Note that the
additive components must have their normalisations frozen for this model to work
properly.

$(FIELDS)

## Example

```julia
model = PowerLaw()
model.K.frozen = true

flux_model = Log10Flux()(model)
```
"""
struct Log10Flux{T} <: AbstractSpectralModel{T, Convolutional}
    E_min::T
    E_max::T
    log10Flux::T
end
function Log10Flux(;
        E_min = FitParam(0.2, lower_limit = 0, frozen = true),
        E_max = FitParam(2.0, lower_limit = 0, frozen = true),
        log10Flux = FitParam(-10.0, lower_limit = -100, upper_limit = 100),
    )
    return Log10Flux(E_min, E_max, log10Flux)
end

function invoke!(flux, energy, model::Log10Flux)
    ilow = clamp(
        _or_else(findfirst(>(model.E_min), energy), length(energy) - 1),
        1,
        lastindex(energy) - 1,
    )
    ihigh = clamp(
        _or_else(findfirst(>(model.E_max), energy), length(energy) - 1) - 1,
        1,
        lastindex(energy) - 1,
    )

    total_e_flux = zero(eltype(flux))

    # low bin straddle
    if ilow > 1
        weight = (energy[ilow]^2 - model.E_min^2) / (energy[ilow] - energy[ilow - 1])
        total_e_flux += flux[ilow - 1] * weight
    end

    for i in ilow:ihigh
        f = flux[i]
        e_low = energy[i]
        e_high = energy[i + 1]

        if e_high > e_low
            total_e_flux += f * (e_high^2 - e_low^2) / (e_high - e_low)
        end
    end

    # high bin straddle
    if ihigh > 1
        weight = (model.E_max^2 - energy[ihigh + 1]^2) / (energy[ihigh + 2] - energy[ihigh + 1])
        total_e_flux += flux[ilow + 1] * weight
    end

    # convert keV to ergs
    total_e_flux = total_e_flux * 0.801096e-9

    flux_exp = 10^model.log10Flux

    if total_e_flux > 0
        @. flux = flux * flux_exp / total_e_flux
    end
    return
end

_or_else(::Nothing, v) = v
_or_else(value, v) = value

export Log10Flux
