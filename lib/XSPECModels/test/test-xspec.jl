using Test
using SpectralFitting
using XSPECModels: @xspecmodel

# just assume that LibXSPEC_jll compiled and shipped okay, and we'll re-wrap a handful of models
# to test that they all work okay with the wrapping macro

@xspecmodel :C_powerlaw struct Test_XS_PowerLaw{T} <: AbstractSpectralModel{T, Additive}
    K::T
    a::T
end
Test_XS_PowerLaw(; K = FitParam(2.0), a = FitParam(1.0)) = Test_XS_PowerLaw(K, a)

model = Test_XS_PowerLaw()

# test invocation
energy = collect(range(0.1, 100.0, 100))
out_flux = invokemodel(energy, model)
# ensure finite
@test all(.!isnan.(out_flux))
@test all(.!isinf.(out_flux))
