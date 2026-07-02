using SpectralFitting, Test


# generate some fake powerlaw data
dummy_data = make_dummy_dataset((E) -> (E^(-3.0)); units = u"counts / (s * keV)")

# test that both julia and xspec implementations can fit simple
model = DummyMultiplicative() * PowerLaw() + PowerLaw() + PowerLaw()

config = SpectralFitting.FittingConfig(FittingProblem(model => dummy_data))

params = get_value.(config.u0)

result = @inferred SpectralFitting.calculate_objective!(config, params)
@test sum(result[1]) ≈ 34.3 atol = 0.01
allocated_bytes = @allocated SpectralFitting.calculate_objective!(config, params)

# TODO: this should be zero
@test allocated_bytes == 48
