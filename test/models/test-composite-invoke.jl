using Test
using SpectralFitting

model = PhotoelectricAbsorption()
domain = collect(range(0.1, 10.0, 100))

# allocating call
out = invokemodel(domain, model)

# Composite models
model = PhotoelectricAbsorption() * PowerLaw()

# allocating call
out = invokemodel(domain, model)
@test sum(out) ≈ 0.501316 atol = 1.0e-3

outputs = allocate_model_output(model, domain)

# in-place
out = invokemodel!(outputs, domain, model)

@test sum(out) ≈ 0.501316 atol = 1.0e-3
