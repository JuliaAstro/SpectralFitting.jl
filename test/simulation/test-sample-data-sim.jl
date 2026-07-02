using Test, SpectralFitting

# path to the data directory
data1 = SpectralFitting.XmmData(joinpath(testdir, "xmm/pn_spec_grp.fits"))

# smoke test
model = GaussianLine() + PowerLaw(a = FitParam(0.2))
sim = simulate(model, data1.response, data1.ancillary; seed = 8, exposure_time = 1.0e1)

# TODO: add a fit. can't do it at the moment as the simulated datasets don't
# support masking the model domain, so we have singular values at high energies
