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

# test write_table_model round-trip
# Create a simple test table model
energy_bins = collect(range(0.1, 10.0, length=21))  # 20 energy bins
param1_grid = [1.0, 2.0, 3.0]
param2_grid = [0.5, 1.0, 1.5]
param_grids = (param1_grid, param2_grid)

n_energy_bins = length(energy_bins) - 1
n_spectra = length(param1_grid) * length(param2_grid)
spectra = zeros(Float64, n_energy_bins, n_spectra)

# Fill with a simple power law model: norm * E^(-index)
let idx = 1
    for p2 in param2_grid
        for p1 in param1_grid
            E_mid = (energy_bins[1:end-1] .+ energy_bins[2:end]) ./ 2
            spectra[:, idx] = p1 .* E_mid .^ (-p2)
            idx += 1
        end
    end
end

# Write the table model
temp_path = tempname() * ".fits"
write_table_model(
    temp_path,
    energy_bins,
    param_grids,
    spectra;
    model_name = "TEST_MODEL",
    param_names = ["norm", "index"],
    param_units = ["", ""],
    model_units = "photons/cm^2/s"
)

@test isfile(temp_path)

# Try to read it back
tmd = TableModelData(Val(2), temp_path; T = Float64)

# Verify energy bins match (approximately, due to float precision)
@test length(tmd.energy_bins) == length(energy_bins)
@test all(isapprox.(tmd.energy_bins, energy_bins, rtol=1e-6))

# Verify parameter grids match
@test length(tmd.params) == 2
@test all(isapprox.(tmd.params[1], param1_grid, rtol=1e-6))
@test all(isapprox.(tmd.params[2], param2_grid, rtol=1e-6))

# Verify spectra match
@test size(tmd.grids) == (length(param1_grid), length(param2_grid))
let idx = 1
    for (j, p2) in enumerate(param2_grid)
        for (i, p1) in enumerate(param1_grid)
            grid_values = tmd.grids[i, j].values
            @test all(isapprox.(grid_values, spectra[:, idx], rtol=1e-6))
            idx += 1
        end
    end
end
