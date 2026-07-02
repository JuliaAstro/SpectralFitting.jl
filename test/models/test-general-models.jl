using Test
using SpectralFitting

# this file tests for specific model constructions
# and anything that produces errors should be added to this
# test set


@ciskip begin
    model = PhotoelectricAbsorption() * (PowerLaw() + BlackBody())

    energy = collect(range(0.2, 8.0, 1000))
    fluxes = zeros(Float64, (length(energy) - 1, objective_cache_count(model)))
    flux = view(fluxes, :, 1)

    # check model invokes okay
    invokemodel!(fluxes, energy, model)
    @test all(.!isnan.(flux))
    @test all(.!isinf.(flux))
end

model = DummyMultiplicative() * DummyMultiplicative() * (DummyAdditive() + DummyAdditive())

checker_all(model)
