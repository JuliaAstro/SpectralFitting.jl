using SpectralFitting, Test, Random

domain = collect(range(0.1, 12, 101))
values = fill(1.0, 100)
high = collect(range(2, 10, 20))

output = zeros(Float64, length(high) - 1)

SpectralFitting.rebin_if_different_domains!(output, high, domain, values)

expected = fill(3.53825, length(high) - 1)
@test output ≈ expected rtol = 1.0e-3


values = (sin.(domain))[1:(end - 1)]

Random.seed!(1)
high = sort(abs.(rand(10))) .* 10 .+ 1.3

output = zeros(Float64, length(high) - 1)

SpectralFitting.rebin_if_different_domains!(output, high, domain, values)

@test output ≈ [
    2.1049690616787626,
    3.1186031150596665,
    2.200909270925761,
    13.951277678105232,
    -19.53623368624383,
    2.5486398574586784,
    5.932901099589066,
    1.2775677281677529,
    0.44147100012051554,
] rtol = 1.0e-3
