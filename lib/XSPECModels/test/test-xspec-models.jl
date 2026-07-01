using Test
using SpectralFitting
using XSPECModels

ALL_XSPEC_MODELS = [
    XS_PowerLaw,
    XS_BlackBody,
    XS_BremsStrahlung,
    XS_Laor,
    XS_DiskLine,
    # XS_KerrDisk, # data file it needs is big
    # XS_KyrLine, # data file it needs is huge
    XS_Gaussian,
    XS_PhotoelectricAbsorption,
    XS_WarmAbsorption,
]

for M in ALL_XSPEC_MODELS
    SpectralFitting.download_model_data(M)
end

ALL_XSPEC_CONVOLUTIONAL = [XS_CalculateFlux]

energy = collect(range(0.1, 100.0, 100))

for M in ALL_XSPEC_MODELS
    m = M()
    # can we invoke okay
    f = invokemodel(energy, m)
    @test all(.!isnan.(f))
    @test all(.!isinf.(f))
end


# have to check convolutional models on a non-zero input flux
for M in ALL_XSPEC_CONVOLUTIONAL
    m = M()
    f = ones(Float64, length(energy) - 1)
    invokemodel!(f, energy, m)
    @test all(.!isnan.(f))
    @test all(.!isinf.(f))
end

# related to additive models should not pass the normalisation constant
# for XS_PowerLaw, a is second argument and the pointer conversion should be offset
# so the shortest test we can make is to set `K` equal and modify `a` to make sure the
# right parameter is passed to the implementation
model1 = XS_PowerLaw(K = FitParam(2.0), a = FitParam(3.0))
model2 = XS_PowerLaw(K = FitParam(2.0), a = FitParam(1.0))
flux1 = invokemodel(energy, model1)
flux2 = invokemodel(energy, model2)
# these should _not_ be the same, despite same normalisation
# so use the sum of residiuals as a metric for difference
@test isapprox(sum(flux1 .- flux2), 86.18438944203572, rtol = 1.0e-4)


# inplace variants
model = XS_PowerLaw()
flux = zeros(Float64, length(energy) - 1)
invokemodel!(flux, energy, model)
@test all(.!isnan.(flux))
@test all(.!isinf.(flux))

# composite
model = XS_PowerLaw() + XS_PowerLaw()
fluxes = hcat(flux, copy(flux))
flux = view(fluxes, :, 1)
invokemodel!(fluxes, energy, model)
@test all(.!isnan.(flux))
@test all(.!isinf.(flux))
