using Test
using SpectralFitting

Base.@kwdef struct ModelTester{T} <: AbstractSpectralModel{T, Additive}
    K::T = FitParam(1.0)
    p1::T = FitParam(3.0)
end

function SpectralFitting.invoke!(out, domain, model::ModelTester)
    @test model.p1 == 5.0
    return
end

function set_to_five(p)
    p.a1.p1 = 5.0
    p.a2.p1 = 5.0
    return
end

model = ParameterPatch(ModelTester() + ModelTester(); patch = set_to_five)

domain = collect(range(0.1, 10.0, 100))
invokemodel(domain, model)

# test the same function can be applied to the model
set_to_five(model)
@test model.a2.p1 == model.a1.p1
@test model.a2.p1.value == 5.0


model = PowerLaw(; K = FitParam(4.0)) + PowerLaw(; K = FitParam(8.0), a = FitParam(8.0))

function relate_normalisations!(p)
    p.a1.K = p.a2.K / 2
    return
end

patched = ParameterPatch(model; patch = relate_normalisations!)
patched.a1.K.frozen = true
patched.a1.K = 0
patched.a2.a.value = 3.0

@test size(SpectralFitting.parameter_vector(patched)) == (4,)
@test SpectralFitting.parameter_names(patched) == (:K, :a, :K, :a)

domain = collect(range(0.1, 10.0, 100))
@inferred invokemodel(domain, patched)

dummy_data = make_dummy_dataset(identity; units = u"counts / (keV * s)")
sim = simulate(model, dummy_data; seed = 42)

prob = FittingProblem(patched => sim)
details(prob)
conf = FittingConfig(prob)

result = fit(prob, LevenbergMarquadt())

apply_patch!(patched)

@test result.u ≈ [2.0004, 7.9989, 8.0001] rtol = 1.0e-3
@test result.stats[1] ≈ 75.639 rtol = 1.0e-3
