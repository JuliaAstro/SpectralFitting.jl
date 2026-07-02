using Test
using SpectralFitting


# multiplicative

model = DummyMultiplicativeTableModel()

@test SpectralFitting.parameter_names(model) == (:a, :b)
@test SpectralFitting.parameter_vector(model) == [model.a, model.b]

# can we invoke the table models alright
energy = collect(range(0.1, 100.0, 100))
out_flux = invokemodel(energy, model)
@test all(out_flux .== 2)

# and compose them
cm = model * model
out_flux = invokemodel(energy, cm)
@test all(out_flux .== 4.0)

# additive

model = DummyAdditiveTableModel()

out_flux = invokemodel(energy, model)
@test all(out_flux .== 3.0)

# composition
cm = model + model
out_flux = invokemodel(energy, cm)
@test all(out_flux .== 6.0)

# compose multiplicative with additive
cm = DummyMultiplicativeTableModel() * DummyAdditiveTableModel()
out_flux = invokemodel(energy, cm)
@test all(out_flux .== 6.0)
