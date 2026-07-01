using Test
using SpectralFitting


# put a couple of delta emission lines together
lines = DeltaLine(; E = FitParam(3.0), K = FitParam(2.0)) + DeltaLine(; E = FitParam(7.0))

# construct the convolutional wrapper
base_model = GaussianLine(; μ = FitParam(1.0), σ = FitParam(0.3))
conv = AsConvolution(base_model)

model = conv(lines)

@test length(SpectralFitting.parameter_vector(model)) == 6
@test SpectralFitting.parameter_names(model) == (:μ, :σ, :K, :E, :K, :E)

domain = collect(range(0.0, 10.0, 150))

output = invokemodel(domain, model)

@test sum(output) ≈ 4.050608829695485 atol = 1.0e-4
@test output[10] ≈ 0.0022170544439135222 atol = 1.0e-4
@test output[40] ≈ 0.058630601782812125 atol = 1.0e-4

# simulate a model spectrum
dummy_data = make_dummy_dataset((E) -> (E^(-3.0)); units = u"counts / (s * keV)")
sim = simulate(model, dummy_data; seed = 42)

prob = FittingProblem(model => dummy_data)

ps, _, bs = SpectralFitting.parameter_vector_symbols_and_bindings(prob)
@test length(ps) == 6

conf = FittingConfig(prob)
output = @inferred SpectralFitting.calculate_objective!(conf, get_value.(conf.u0))

model.c1.μ.frozen = true
model.a1.K.frozen = true
model.a2.K.frozen = true
model.a1.E.frozen = true
model.a2.E.frozen = true

# change the width
model.c1.σ.value = 0.1
model

prob = FittingProblem(model => sim)
conf = FittingConfig(prob)

@test length(prob.lookup) == 6

@inferred SpectralFitting.calculate_objective!(conf, [2.0])

result = fit(prob, LevenbergMarquadt(), verbose = true)

@test sum(result.stats) ≈ 76.71272868245076 atol = 1.0e-3
@test result.u[1] ≈ 0.3 atol = 1.0e-2

# put a couple of delta emission lines together
lines =
    DeltaLine(; E = FitParam(3.0), K = FitParam(2.0), width = 0.1) +
    DeltaLine(; E = FitParam(7.0))
model = conv(lines)

sim = simulate(model, dummy_data; seed = 42)

# now see if we can fit the delta line
model.c1.μ.frozen = true
model.a1.K.frozen = true
model.a2.K.frozen = true
model.a1.E.frozen = true
model.a2.E.frozen = true
model.c1.σ.frozen = true

model.a1.E.frozen = false
model.a1.E.value = 2.0

model
begin
    prob = FittingProblem(model => sim)
    result = fit(prob, LevenbergMarquadt(); verbose = true)
end
@test sum(result.stats) ≈ 76.66970981760741 atol = 1.0e-3
@test result.u[1] ≈ 3.0 atol = 1.0e-2
