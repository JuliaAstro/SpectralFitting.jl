using Test, SpectralFitting



dummy_data = make_dummy_dataset((E) -> (E^(-3.0)); units = u"counts / (s * keV)")
model = PowerLaw()
prob = FittingProblem(model, dummy_data)

# test inference
@inferred SpectralFitting._unpack_config(prob);

result = fit(prob, LevenbergMarquadt())

@test result[1].u ≈ [12.066, 3.0810] atol = 1e-3

SpectralFitting.calculate_objective!(result[1], result[1].u)

# test API
SpectralFitting.get_objective(result[1])
SpectralFitting.get_objective_variance(result[1])

@test measure(ChiSquared(), result) ≈ sum(result.stats) rtol = 1e-3
@test measure(ChiSquared(), result, [11.0, 2.0]) ≈ 281423.5 rtol = 1e-3
# can measure different parameters
@test measure(ChiSquared(), result[1], [11.0, 2.0]) ≈ 281423.5 rtol = 1e-3

# check these don't error
@test calculate_objective!(result, result.u)[1] ≈
      calculate_objective!(result[1], result[1].u)
calculate_objective!(result, [11.0, 2.0])

# mismatch
@test_throws "" calculate_objective!(result, [11.0])
@test_throws "" calculate_objective!(result[1], [11.0])

# test reduced statistic
# dof = n_bins - n_free_params = 100 - 2 = 98
# reduced_statistic should be χ² (stats) / dof (98)
@test dof(result[1]) == 98
@test reduced_statistic(result[1]) ≈ result[1].stats / 98 rtol = 1e-10

# multi-dataset with one tied parameter (K)
dummy_data2 = make_dummy_dataset((E) -> (E^(-3.0)); units = u"counts / (s * keV)")
prob_tied = FittingProblem(PowerLaw() => dummy_data, PowerLaw() => dummy_data2)
bindall!(prob_tied, :K)
result_tied = fit(prob_tied, LevenbergMarquadt())

# K is shared, so there are only 3 free params for all
@test count(result_tied.config.parameter_cache.free_mask) == 3

# dof for each slice shows that K still counts as free for each dataset
@test dof(result_tied[1]) == 98
@test dof(result_tied[2]) == 98

# use the free_mask directly instead of dof function to check total free params across all datasets 
# should not double count the tied parameter, so total free params should be 3, not 4
total_n_bins =
    length(result_tied.config.data_cache[1].objective) +
    length(result_tied.config.data_cache[2].objective)
@test total_n_bins - count(result_tied.config.parameter_cache.free_mask) == 197