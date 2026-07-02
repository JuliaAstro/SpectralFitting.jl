using SpectralFitting, Test


data = make_dummy_dataset(collect(range(1.0, 0.0, 10)), collect(range(0, 15.0, 11)))

# this really needs to be compile time known
@inferred preferred_units(data, ChiSquared())
@inferred preferred_units(typeof(data), ChiSquared())

set_units!(data, u"counts / (s * keV)")

make_objective(SpectralFitting.ContiguouslyBinned(preferred_units(data, Cash())), data)
