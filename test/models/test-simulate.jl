using Test
using SpectralFitting


dummy_data = make_dummy_dataset((E) -> (E^(-3.0)); units = u"counts / (s * keV)")
model = PowerLaw()

sim = simulate(model, dummy_data; seed = 42)

result = fit(FittingProblem(model => sim), LevenbergMarquadt())

@test result.u ≈ [1.0001, 2.0003] atol = 1.0e-4
@test result.stats[1] ≈ 76.923 atol = 1.0e-3

simulate!(result.config, [3.0, 0.1])

result2 = fit(result.config, LevenbergMarquadt())

@test result2.u ≈ [2.99, 0.1] atol = 1.0e-2
@test result2.stats[1] ≈ 95.519 atol = 1.0e-3

# since goodness has wrappers around simulate, we'll do a quick test here too

pcent, trials = goodness(result2[1]; seed = 42, N = 100)

@test pcent ≈ 41.0 atol = 1.0e-2
