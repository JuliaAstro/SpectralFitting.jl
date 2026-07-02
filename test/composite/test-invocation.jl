using Test
using SpectralFitting


include("../fuzz.jl")

model = DummyAdditive()

# test single invocation
energy = collect(range(0.1, 100.0, 100))
flux = zeros(Float64, length(energy) - 1)

out_flux = invokemodel!(flux, energy, model)
@test all(out_flux .== 6)

# check type stability
@inferred invokemodel!(flux, energy, model)

# make sure normalisation is handled correctly
model2 = DummyAdditive(K = FitParam(2.0))

out_flux = invokemodel!(flux, energy, model2)
@test all(out_flux .== 12)

# just check that it has no errors
checker_all(model) = invokemodel(energy, model)
_ = checker_all(model + model)

# some other special cases
_ = checker_all(DummyMultiplicative() * DummyAdditive())

# added models
out_flux = invokemodel(energy, model + model)
@test all(out_flux .== 12)
@inferred invokemodel(energy, model + model)

# big added models
out_flux = @inferred invokemodel(
    energy,
    model +
        model +
        model +
        model +
        model +
        model +
        model +
        model +
        model +
        model +
        model +
        model +
        model +
        model +
        model +
        model,
)
@test all(out_flux .== 96)

# test multiplicative with additive
cm = DummyMultiplicative() * DummyAdditive()
out_flux = invokemodel(energy, cm)
@test all(out_flux .== 30)

# test inplace version
model = DummyAdditive()
invokemodel!(out_flux, energy, model)
@test all(out_flux .== 6)

# should be unchanged if called multiple times
@inferred invokemodel!(out_flux, energy, model)
@inferred invokemodel!(out_flux, energy, model)
@test all(out_flux .== 6)

# composite with many frozen
cm = @inferred DummyMultiplicative() * (DummyAdditive() + DummyAdditiveWithManyFrozen())
flux = invokemodel(energy, cm)
@test all(flux .== 160.0)

# inplace variant
fluxes = allocate_model_output(cm, energy)
flux = view(fluxes, :, 1)
invokemodel!(fluxes, energy, cm)
@test all(flux .== 160.0)

# edge case that previously was failing
model = DummyMultiplicative() * DummyMultiplicative() * (DummyAdditive() + DummyAdditive())
flux = invokemodel(energy, model)
@test all(flux .== 300)

fluxes = zeros(Float64, (length(energy) - 1, objective_cache_count(model)))
flux = view(fluxes, :, 1)
invokemodel!(fluxes, energy, model)
@test all(flux .== 300)

# assert error is raised when too few flux passed
@test_throws "" invokemodel!(flux, energy, model)
