using Test
using SpectralFitting


# generate some fake powerlaw data with three components
dummy_data = make_dummy_dataset(
    (E) -> (E^(-0.1) + E^(-3.0) + E^(-1.0));
    units = u"counts / (s * keV)",
)

# model with two components
model = PowerLaw() + PowerLaw()
prob = FittingProblem(
    FittableMultiModel(model, model),
    FittableMultiDataset(dummy_data, dummy_data),
)
# test apis all work
prob = FittingProblem(model => dummy_data, model => dummy_data)

result = fit(prob, LevenbergMarquadt())

res1 = result[1]
res2 = result[2]

# both models should fit more or less the same
@test res1.u ≈ res2.u atol = 1.0e-2
@test res1.stats ≈ res2.stats atol = 1.0e-2

# now change second models data to ensure the normalisations fit independelty
dummy_data2 = deepcopy(dummy_data)
dummy_data2.spectrum.data .*= 3.0

prob = FittingProblem(
    FittableMultiModel(model, model),
    FittableMultiDataset(dummy_data, dummy_data2),
)

result = fit(prob, LevenbergMarquadt())
# photon indices should still be the same
r1 = result[1]
r2 = result[2]
@test r1.u[2] + r1.u[4] ≈ r2.u[2] + r2.u[4]
# but norms should be different
@test !(r1.u[1] + r1.u[3] ≈ r2.u[1] + r2.u[3])

# now we bind one of the normalisations
bindall!(prob, (:a1, :K))

result = fit(prob, LevenbergMarquadt())
# ensure bound parameter is identical
r1 = result[1]
r2 = result[2]
@test r1.u[1] == r2.u[1]
@test !(r1.u[3] ≈ r2.u[3])

# now the second model will be different
model2 = PowerLaw() + PowerLaw() + PowerLaw()

prob = FittingProblem(
    FittableMultiModel(model, model2),
    FittableMultiDataset(dummy_data, dummy_data2),
)

# ensure we can fit this fine
result = fit(prob, LevenbergMarquadt())

@test result[1].stats ≈ 2.86 atol = 1.0e-2
@test result[2].stats ≈ 0.079 atol = 1.0e-3

# and we can still bind
prob = FittingProblem(
    FittableMultiModel(model, model2),
    FittableMultiDataset(dummy_data, dummy_data2),
)
bindall!(prob, (:a1, :K))

result = fit(prob, LevenbergMarquadt())

@test result[1].stats ≈ 2.8627679111508226 atol = 0.1
@test result[2].stats ≈ 0.2460801839134885 atol = 0.1
