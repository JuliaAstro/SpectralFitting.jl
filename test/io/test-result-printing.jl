using Test
using SpectralFitting

function showstring(item)
    buffer = IOBuffer()
    Base.show(IOContext(buffer, :color => false), MIME"text/plain"(), item)
    strip(String(take!(buffer)))
end

# generate some fake powerlaw data
model = PowerLaw()
dummy_data = make_dummy_dataset((E) -> (E^(-3.0)); units = u"counts / (s * keV)")
result = fit(FittingProblem(model, dummy_data), LevenbergMarquadt())

# test result printing is vertical with labelled names and has correct values
slice_string = showstring(result[1])

# Expected output (formatted):
#
# ┌ FitResultSlice:
# │ Model: PowerLaw
# │  . Name    : u         Δu
# │  . m.K     : 12.066    0.24374
# │  . m.a     : 3.0810    0.012679
# │
# │  . χ²      : 14.963
# │  . χᵥ²     : 0.15268 (dof=98)
# └
expected_slice = "┌ FitResultSlice:\n│ Model: PowerLaw\n│  . Name    : u         Δu        \n│  . m.K     : 12.066    0.24374   \n│  . m.a     : 3.0810    0.012679  \n│ \n│  . χ²      : 14.963\n│  . χᵥ²     : 0.15268 (dof=98)\n└"
@test slice_string == expected_slice

# Two-dataset with K bound across models: global dof should be 3 and use free_mask count, not sum of
# per-slice counts (which would double-count the shared K and give dof 196 instead of 197)
dummy_data2 = make_dummy_dataset((E) -> (E^(-3.0)); units = u"counts / (s * keV)")
prob_bound = FittingProblem(PowerLaw() => dummy_data, PowerLaw() => dummy_data2)
bindall!(prob_bound, :K)
result_bound = fit(prob_bound, LevenbergMarquadt())
result_bound_string = showstring(result_bound)

# Expected output (formatted):
#
# ┌ FitResult:
# │  Model: PowerLaw
# │   . Name    : u         Δu
# │   . m.K     : 12.066    0.17235
# │   . m.a     : 3.0810    0.010006
# │
# │   . χ²      : 14.963
# │   . χᵥ²     : 0.15268 (dof=98)
# │  Model: PowerLaw
# │   . Name    : u         Δu
# │   . m.K     : 12.066    0.17235
# │   . m.a     : 3.0810    0.010006
# │
# │   . χ²      : 14.963
# │   . χᵥ²     : 0.15268 (dof=98)
# └ Σχ² = 29.925, χᵥ² = 0.15190 (dof=197)
#
# total dof = 200 total bins - 3 free params (K shared, a in first dataset, a in second dataset) = 197
expected_result_bound = "┌ FitResult:\n│  Model: PowerLaw\n│   . Name    : u         Δu        \n│   . m.K     : 12.066    0.17235   \n│   . m.a     : 3.0810    0.010006  \n│  \n│   . χ²      : 14.963\n│   . χᵥ²     : 0.15268 (dof=98)\n│  Model: PowerLaw\n│   . Name    : u         Δu        \n│   . m.K     : 12.066    0.17235   \n│   . m.a     : 3.0810    0.010006  \n│  \n│   . χ²      : 14.963\n│   . χᵥ²     : 0.15268 (dof=98)\n└ Σχ² = 29.925, χᵥ² = 0.15190 (dof=197)"
@test result_bound_string == expected_result_bound
