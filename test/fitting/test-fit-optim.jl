using Test
using SpectralFitting
using OptimizationOptimJL


# generate some fake powerlaw data with three components
dummy_data = make_dummy_dataset(
    (E) -> (E^(-0.1) + E^(-3.0) + E^(-1.0));
    units = u"counts / (s * keV)",
)

# model with two components
model = PowerLaw(K = FitParam(10.0)) + PowerLaw(K = FitParam(10.0))
prob = FittingProblem(
    FittableMultiModel(model, model),
    FittableMultiDataset(dummy_data, dummy_data),
)

# check inferred types for multi models
@inferred SpectralFitting._unpack_config(prob)

result = fit(prob, NelderMead(), autodiff = SpectralFitting.Optimization.SciMLBase.NoAD())

# both models should fit more or less the same
@test sum(result[1].u) ≈ sum(result[2].u) atol = 1.0e-2
@test result[1].stats ≈ result[2].stats atol = 1.0e-2
@test sum(result.stats) ≈ 5.7238519700611725 atol = 1.0e-2

# now with the AD backend
result = fit(prob, BFGS())

# both models should fit more or less the same
@test sum(result[1].u) ≈ sum(result[2].u) atol = 1.0e-2
@test result[1].stats ≈ result[2].stats atol = 1.0e-2
@test sum(result.stats) ≈ 5.7238519700611725 atol = 1.0e-2

# do a single model

prob = FittingProblem(model => dummy_data)

result = fit(prob, BFGS())

# both models should fit more or less the same
@test sort(result.u) ≈ sort([17.254, 0.285, 14.97, 3.0974]) atol = 1.0e-2
@test sum(result.stats) ≈ 2.8619 atol = 1.0e-2

# now with different statistic

prob = FittingProblem(model => dummy_data)

result = fit(prob, BFGS(); stat = Cash())

# both models should fit more or less the same
@test sort(result.u) ≈ sort([19.047, 0.339, 12.207, 3.289]) atol = 1.0e-2
@test sum(result.stats) ≈ 0.305 atol = 1.0e-2
