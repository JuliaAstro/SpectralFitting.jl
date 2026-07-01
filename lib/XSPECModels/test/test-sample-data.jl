using Test
using SpectralFitting
using XSPECModels

# might have to download the model data first
SpectralFitting.download_model_data(XS_Laor)
SpectralFitting.download_model_data(PhotoelectricAbsorption)

function prepare_data!(data, low, high)
    drop_bad_channels!(data)
    regroup!(data)
    normalize!(data)
    mask_energies!(data, low, high)
    return data
end

# path to the data directory
data1 = SpectralFitting.XmmData(joinpath(testdir, "xmm/pn_spec_grp.fits"))
prepare_data!(data1, 0.8, 10.0)

model = PhotoelectricAbsorption() * XS_PowerLaw() + XS_Laor()

# construct the model and data problem
prob = FittingProblem(model, data1)

result = fit(prob, LevenbergMarquadt())

result[1].u

# todo: check the chi2 calculations
@test sum(result.stats) ≈ 484.5 atol = 0.1

@test result.u ≈ [0.1376, 0.010446, 1.8422, 8.9215e-5, 6.5673] rtol = 1.0e-3
xspec_u = [0.140182, 1.048e-2, 1.84862, 9.25e-5, 6.559]
@test result.u ≈ xspec_u rtol = 1.0e-2

data2 = NuStarData(joinpath(testdir, "nustar/nu60001047002A01_sr_grp_simple.pha"))
prepare_data!(data2, 2.0, 14.0)

prob = FittingProblem(model, data2)
result = fit(prob, LevenbergMarquadt())

@test sum(result.stats) ≈ 191.0 atol = 0.1
@test result.u ≈ [
    0.33187746370113763,
    0.01989636801877023,
    1.9972614257180774,
    0.0002059507411466161,
    6.499205973061382,
] rtol = 1.0e-2

# test joint
prob = FittingProblem(model => data1, model => data2)
result = fit(prob, LevenbergMarquadt())

@test result[1].stats ≈ 484.546 atol = 0.1
@test result[2].stats ≈ 190.999 atol = 0.1

# with binding
bindall!(prob, (:a2, :K))

result = fit(prob, LevenbergMarquadt())

@test result[1].stats ≈ 485.0180163473964 atol = 0.1
@test result[2].stats ≈ 202.1704312745641 atol = 0.1

# todo: with background subtraction
data1_nobkg = deepcopy(data1)
set_units!(data1_nobkg, u"counts")
subtract_background!(data1_nobkg)
set_units!(data1_nobkg, u"counts / (s * keV)")

prob = FittingProblem(model, data1_nobkg)
result = fit(prob, LevenbergMarquadt())
# these have been checked and are the same as XSPEC
@test sum(result.stats) ≈ 496.0927976691264 atol = 0.1
@test result.u ≈ [0.1395103, 0.01047537, 1.848017, 9.18112e-5, 6.55961] rtol = 1.0e-3

# try different domain entirely
set_domain!(data1, collect(range(0.01, 20.0, 1000)))

model = PhotoelectricAbsorption() * XS_PowerLaw() + XS_Laor()

# construct the model and data problem
prob = FittingProblem(model, data1)

result = fit(prob, LevenbergMarquadt())
# these have been checked and are the same as XSPEC
@test sum(result.stats) ≈ 466.6477135605372 atol = 0.1
@test result.u ≈ [0.13429082, 0.010527296, 1.8470323, 9.076584e-5, 6.574072] rtol = 1.0e-3
