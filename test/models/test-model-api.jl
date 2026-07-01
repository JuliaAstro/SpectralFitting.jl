using Test, SpectralFitting


# single

model = DummyAdditive()
domain = collect(range(0.1, 10, 200))
output = @inferred allocate_model_output(model, domain)
@test size(output) == (199, 1)

invokemodel!(output, domain, model)
@test all(i -> isapprox(i, 6.0), output)

output = invokemodel(domain, model)
@test all(i -> isapprox(i, 6.0), output)

invokemodel!(output, domain, model)
allocated = @allocated invokemodel!(output, domain, model)
@test allocated <= 209

params = get_value.(SpectralFitting.parameter_vector(model))
invokemodel!(output, domain, model, params)
allocated = @allocated invokemodel!(output, domain, model, params)
@test allocated <= 48

@test modelkind(model) == Additive()

# composite

model = DummyAdditive() + DummyAdditive()
output = @inferred allocate_model_output(model, domain)
@test size(output) == (199, 2)

invokemodel!(output, domain, model)
@test all(i -> isapprox(i, 12.0), output[:, 1])

invokemodel!(output, domain, model)

params = get_value.(SpectralFitting.parameter_vector(model))
invokemodel!(output, domain, model, params)
allocated = @allocated invokemodel!(output, domain, model, params)
@test allocated <= 48

output = invokemodel(domain, model)
@test all(i -> isapprox(i, 12.0), output)

@test modelkind(model) == Additive()

# parameter things

@test SpectralFitting.parameter_count(PhotoelectricAbsorption()) == 1
@test SpectralFitting.parameter_count(PowerLaw()) == 2
@test SpectralFitting.parameter_count(PowerLaw() + PowerLaw()) == 4

model = DummyMultiplicative() * PowerLaw(a = FitParam(4.0, frozen = true))
cache = @inferred SpectralFitting.make_parameter_cache(model)
SpectralFitting.update_free_parameters!(cache, [2.0, 0.5])
@test cache.parameters == [2.0, 5.0, 0.5, 4.0]

outputs = allocate_model_output(model, domain)
output = invokemodel!(outputs, domain, model, cache)
@test sum(output) ≈ 1666.665 atol = 1.0e-4

model = DummyAdditive()

params_tuple = @inferred SpectralFitting._unpack_as_tuple(model)
expected = (FitParam(1.0), FitParam(1.0), FitParam(5.0))
@test isapprox.(params_tuple, expected) |> all
