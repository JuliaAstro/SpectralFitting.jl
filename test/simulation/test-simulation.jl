using Test, SpectralFitting


dummy_data = make_dummy_dataset((E) -> (E^(-3.0)); units = u"counts / (s * keV)")
model = PowerLaw()

sim_data = simulate(model, dummy_data; seed = 2, exposure_time = 1.0e5)

@test sum(sim_data.data) ≈ 4.897 atol = 1.0e-2

prob = FittingProblem(model => sim_data)
result = fit(prob, LevenbergMarquadt())

@test result.u ≈ [1.0, 2.0] atol = 1.0e-2

# TODO: test the simulation api
sim_data2 = simulate(model, dummy_data.response; seed = 2)
@test sum(sim_data2.data) ≈ 4.897 atol = 1.0e-2
