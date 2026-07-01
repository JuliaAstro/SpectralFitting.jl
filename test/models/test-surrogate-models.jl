using Test
using SpectralFitting
import Random

lower_bounds = (1.0e-3,)
upper_bounds = (3.0,)

energy = collect(range(0.1, 20.0, 100))
model = PhotoelectricAbsorption()

flux = similar(energy)[1:(end - 1)]

surrogate = make_surrogate_harness(
    (x, y) -> SpectralFitting.RadialBasis(x, y, lower_bounds, upper_bounds),
    energy,
    model,
    lower_bounds,
    upper_bounds,
)

# number of points the surrogate has been trained on
length(surrogate.surrogate.x)

nh = 2.0
model.ηH.value = nh

f = invokemodel(energy, model)

f̂ = surrogate.surrogate([nh])

@test sum(f .- f̂) ≈ 0.0 atol = 1.0e-3

optimize_accuracy!(surrogate; maxiters = 50)

f̂ = surrogate.surrogate([nh])

@test sum(f .- f̂) ≈ -0.000362214 atol = 1.0e-4

sm = @inferred make_model(surrogate)

# smoke test these
invokemodel(energy, sm)

new_model = sm * PowerLaw()
invokemodel(energy, new_model)
