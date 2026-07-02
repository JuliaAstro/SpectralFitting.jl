using Test, SpectralFitting


# composite power law
dummy_data = make_dummy_dataset((E) -> (E^(-3.0) + E^(-1)); units = u"counts / (s * keV)")

# julia implementation
model = PowerLaw() + PowerLaw(a = FitParam(1.0))

model.a1.K = 12
model.a1.K.frozen = true

model

prob = FittingProblem(model, dummy_data)
result = fit(prob, LevenbergMarquadt())

@test result.u ≈ [3.2462951946372813, 10.404292824678521, 1.0] rtol = 1.0e-2

update_model!(model, result)

# make sure the correct parameters have been updated
@test model.a1.K.value == 12.0
@test model.a1.a.value ≈ 3.246 rtol = 1.0e-2
@test model.a2.K.value ≈ 10.404 rtol = 1.0e-2
