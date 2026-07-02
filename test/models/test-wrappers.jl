using Test
using SpectralFitting


dataset = make_dummy_dataset(i -> i^-2)

struct TestModelWrapper{M, T, K} <: SpectralFitting.AbstractModelWrapper{M, T, K}
    model::M
end

TestModelWrapper(model::AbstractSpectralModel{T, K}) where {T, K} =
    TestModelWrapper{typeof(model), T, K}(model)

# single wrapper
model = TestModelWrapper(PowerLaw())
domain = collect(range(0.1, 10.0, 100))
@inferred invokemodel(domain, model)

# composite wrapper
model = TestModelWrapper(PowerLaw() + PowerLaw())
domain = collect(range(0.1, 10.0, 100))

output = @inferred invokemodel(domain, model)
@test size(allocate_model_output(model, domain)) == (99, 2)
@test size(SpectralFitting.parameter_vector(model)) == (4,)

backing_output = invokemodel(domain, SpectralFitting.backing_model(model))

@test output ≈ backing_output rtol = 0.01

@inferred SpectralFitting.parameter_names(typeof(model))

prob = FittingProblem(model => dataset)
prob_check = FittingProblem(SpectralFitting.backing_model(model) => dataset)

@test length(prob.lookup) == length(prob_check.lookup)

conf = @inferred FittingConfig(prob)
conf_check = @inferred FittingConfig(prob)

a1 = SpectralFitting.calculate_objective!(conf, [1.0, 4.0, 3.0, 2.0]) |> first
a2 = SpectralFitting.calculate_objective!(conf_check, [1.0, 4.0, 3.0, 2.0]) |> first

@test a1 ≈ a2 rtol = 1.0e-4
@test SpectralFitting._get_parameters(conf.parameter_cache, 1) ≈
    SpectralFitting._get_parameters(conf_check.parameter_cache, 1) atol = 1.0e-3

# some smoke tests
@inferred SpectralFitting.calculate_objective!(conf, [1.0, 1.0, 1.0, 1.0])

result = fit(prob, LevenbergMarquadt())
result_check = fit(prob_check, LevenbergMarquadt())

@test result.u ≈ result_check.u rtol = 1.0e-1
@test result.stats ≈ result_check.stats rtol = 1.0e-1

# do frozen parameters work also?
model.a1.a.frozen = true
prob = FittingProblem(model => dataset)
result = fit(prob, LevenbergMarquadt())

@test result.u ≈ [5.114, 0.51301, 3.0971] atol = 1.0e-2

# composite wrapper but as part of another composite model
model = TestModelWrapper(PowerLaw() + PowerLaw()) + PowerLaw()
domain = collect(range(0.1, 10.0, 100))
@inferred invokemodel(domain, model)
