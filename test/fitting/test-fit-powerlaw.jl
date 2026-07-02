using Test
using SpectralFitting


# generate some fake powerlaw data
dummy_data = make_dummy_dataset((E) -> (E^(-3.0)); units = u"counts / (s * keV)")

# test that both julia and xspec implementations can fit simple
model = PowerLaw()
prob = FittingProblem(model, dummy_data)

# test inference
@inferred SpectralFitting._unpack_config(prob);

result = fit(prob, LevenbergMarquadt())

@test result.u ≈ [12.06629478087094, 3.080992500319396] atol = 1.0e-4

# composite power law
dummy_data = make_dummy_dataset((E) -> (E^(-3.0) + E^(-1)); units = u"counts / (s * keV)")

# julia implementation
model = PowerLaw() + PowerLaw(a = FitParam(1.0))
prob = FittingProblem(model, dummy_data)
result = fit(prob, LevenbergMarquadt())
# for these the order keeps going wrong on the CI, so we'll just check the sum
@test sort(result.u) ≈
    sort([10.374890467986033, 1.0, 12.277664914150247, 3.225587384436137]) atol = 1.0e-1
