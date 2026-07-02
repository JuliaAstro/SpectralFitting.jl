using Test
using SpectralFitting


# generate some fake powerlaw data
dummy_data = make_dummy_dataset((E) -> (E^(-3.0)); units = u"counts / (s * keV)")

# test that both julia and xspec implementations can fit simple
model = PowerLaw()
prob = FittingProblem(model, dummy_data)

result = fit(prob, LevenbergMarquadt())

append!(prob.data.extension[1].low, logrange(1.0e-4, 0.9, 100))
append!(prob.data.extension[1].high, logrange(10.1, 20.0, 100))

conf = FittingConfig(prob)

@test length(conf.data_cache[1].model_domain) == 301
@test extrema(conf.data_cache[1].model_domain) == (1.0e-4, 20.0)

ext_result = fit(prob, LevenbergMarquadt())

@test ext_result.u ≈ result.u rtol = 1.0e-4
