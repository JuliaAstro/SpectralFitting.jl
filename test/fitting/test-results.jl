using Test, SpectralFitting


dummy_data = make_dummy_dataset((E) -> (E^(-3.0)); units = u"counts / (s * keV)")
model = PowerLaw()
prob = FittingProblem(model, dummy_data)

# test inference
@inferred SpectralFitting._unpack_config(prob);

result = fit(prob, LevenbergMarquadt())

@test result[1].u ≈ [12.066, 3.081] atol = 1.0e-3

SpectralFitting.calculate_objective!(result[1], result[1].u)

# test API
SpectralFitting.get_objective(result[1])
SpectralFitting.get_objective_variance(result[1])

@test measure(ChiSquared(), result) ≈ sum(result.stats) rtol = 1.0e-3
@test measure(ChiSquared(), result, [11.0, 2.0]) ≈ 281423.5 rtol = 1.0e-3
# can measure different parameters
@test measure(ChiSquared(), result[1], [11.0, 2.0]) ≈ 281423.5 rtol = 1.0e-3

# check these don't error
@test calculate_objective!(result, result.u)[1] ≈
    calculate_objective!(result[1], result[1].u)
calculate_objective!(result, [11.0, 2.0])

# mismatch
@test_throws "" calculate_objective!(result, [11.0])
@test_throws "" calculate_objective!(result[1], [11.0])
