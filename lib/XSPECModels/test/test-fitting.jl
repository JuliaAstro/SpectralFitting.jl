using Test
using SpectralFitting
using XSPECModels

include("../../../test/dummies.jl")

# generate some fake powerlaw data
dummy_data = make_dummy_dataset(E -> E^(-3.0); units = u"counts / (s * keV)")

model = XS_PowerLaw()
prob = FittingProblem(model, dummy_data)
result = fit(prob, LevenbergMarquadt())
@test result.u ≈ [12.06629478087094, 3.080992500319396] atol = 1.0e-4

# composite power law
dummy_data = make_dummy_dataset(E -> (E^(-3.0) + E^(-1)); units = u"counts / (s * keV)")

# xpsec implementation
model = XS_PowerLaw() + XS_PowerLaw(a = FitParam(1.0))
prob = FittingProblem(model, dummy_data)
result = fit(prob, LevenbergMarquadt())
@test sort(result.u) ≈
    sort([11.152658761725311, 1.0372890868418811, 11.219122689264507, 3.310562124158435]) atol =
    1.0e-1

# mixed implementation
model = XS_PowerLaw() + PowerLaw(a = FitParam(1.0))
prob = FittingProblem(model, dummy_data)
result = fit(prob, LevenbergMarquadt())
@test sort(result.u) ≈
    sort([11.219123015018498, 3.310562096287135, 11.152658542859502, 1.0372890769565792]) atol =
    1.0e-1
