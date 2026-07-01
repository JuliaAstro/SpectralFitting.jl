"""
    XS_CalculateFlux(E_min, E_max, lg10Flux)

$(FIELDS)
"""
@xspecmodel :C_cflux struct XS_CalculateFlux{T} <: AbstractSpectralModel{T, Convolutional}
    "Minimum energy."
    E_min::T
    "Maximum energy."
    E_max::T
    "log (base 10) flux in erg / cm^2 / s"
    log10Flux::T
end
function XS_CalculateFlux(;
        E_min = FitParam(0.2, frozen = true),
        E_max = FitParam(2.0, frozen = true),
        log10Flux = FitParam(-10.0, lower_limit = -100, upper_limit = 100),
    )
    return XS_CalculateFlux(E_min, E_max, log10Flux)
end

"""
    XS_Kerrconv()

$(FIELDS)
"""
@xspecmodel :C_kerrconv struct XS_Kerrconv{T} <: AbstractSpectralModel{T, Convolutional}
    "The emissivity index for the inner disk."
    Index1::T
    "The emissivity index for the outer disk."
    Index2::T
    "The break radius separating the inner and outer portions of the disk, in gravitational radii."
    r_br_g::T
    "The dimensionless spin parameter of the black hole."
    a::T
    "The disk inclination angle, in degrees. A face-on disk has Incl=0."
    Incl::T
    "The inner radius of the disk, in units of the radius of marginal stability."
    Rin_ms::T
    "The outer radius of the disk, in units of the radius of marginal stability."
    Rout_ms::T
end
function XS_Kerrconv(;
        Index1 = FitParam(3.0, frozen = true, lower_limit = -10.0, upper_limit = 10.0),
        Index2 = FitParam(3.0, frozen = true, lower_limit = -10.0, upper_limit = 10.0),
        r_br_g = FitParam(6.0, frozen = true, lower_limit = 1.0, upper_limit = 400.0),
        a = FitParam(0.0, frozen = false, lower_limit = 0.0, upper_limit = 0.998),
        Incl = FitParam(30.0, frozen = false, lower_limit = 0.0, upper_limit = 90.0),
        Rin_ms = FitParam(1.0, frozen = true, lower_limit = 1.0, upper_limit = 400.0),
        Rout_ms = FitParam(400.0, frozen = true, lower_limit = 1.0, upper_limit = 400.0),
    )
    return XS_Kerrconv(Index1, Index2, r_br_g, a, Incl, Rin_ms, Rout_ms)
end

export XS_CalculateFlux, XS_Kerrconv
