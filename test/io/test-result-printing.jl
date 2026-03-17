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
# └
expected_slice = "┌ FitResultSlice:\n│ Model: PowerLaw\n│  . Name    : u         Δu        \n│  . m.K     : 12.066    0.24374   \n│  . m.a     : 3.0810    0.012679  \n│ \n│  . χ²      : 14.963\n└"
@test slice_string == expected_slice
