using SpectralFitting
using Test
using Random

function _make_dataset(K, a, energy_bins; seed=1)
    Random.seed!(seed)
    m    = PowerLaw(K=FitParam(K), a=FitParam(a))
    flux = collect(invokemodel(energy_bins, m))
    var  = abs.(flux) .* 0.001
    InjectiveData(energy_bins[1:end-1], flux; codomain_variance=var)
end

const EBINS  = collect(range(0.3, 10.0, length=51))
const TRUE_K = [1.0, 0.5, 2.0]
const TRUE_a = 1.7

@testset "domain_union" begin
    d1 = _make_dataset(1.0, TRUE_a, EBINS; seed=1)
    d2 = _make_dataset(0.5, TRUE_a, EBINS; seed=2)
    d3 = _make_dataset(2.0, TRUE_a, EBINS; seed=3)

    u = domain_union(d1, d2, d3)
    @test length(u) == 51
    @test first(u) ≈ 0.3  atol=1e-6
    @test last(u)  ≈ 10.0 atol=1e-6
    @test issorted(u)
    @test domains_are_compatible(d1, d2, d3)
    @test_throws ArgumentError domain_union(d1)
end

@testset "extract_normalisation_index" begin
    m = PowerLaw()
    @test extract_normalisation_index(m) == 1
    @test has_normalisation(m) == true
    @test SpectralFitting.parameter_names(m)[1] == :K
    @test extract_normalisation_index(m) !== nothing
end

@testset "multi_dataset_problem structure" begin
    datasets = [_make_dataset(TRUE_K[i], TRUE_a, EBINS; seed=i)
                for i in 1:3]
    base = PowerLaw(
        K = FitParam(1.0, lower_limit=0.0, upper_limit=100.0),
        a = FitParam(1.5, lower_limit=0.0, upper_limit=5.0),
    )
    prob = multi_dataset_problem(base, datasets...)

    @test SpectralFitting.model_count(prob) == 3
    @test SpectralFitting.data_count(prob)  == 3

    # `a` is shared — bindings dict should have entries
    @test !isempty(prob.bindings)

    # Error on single dataset
    @test_throws ArgumentError multi_dataset_problem(base, datasets[1])
end

@testset "multi_dataset_problem fitting" begin
    datasets = [_make_dataset(TRUE_K[i], TRUE_a, EBINS; seed=i)
                for i in 1:3]
    base = PowerLaw(
        K = FitParam(1.0, lower_limit=0.0, upper_limit=100.0),
        a = FitParam(1.5, lower_limit=0.0, upper_limit=5.0),
    )
    prob   = multi_dataset_problem(base, datasets...)
    result = fit(prob, LevenbergMarquadt())

    for i in 1:3
        @test isapprox(result[i].u[1], TRUE_K[i], rtol=0.05)
    end
    @test isapprox(result[1].u[2], TRUE_a, rtol=0.05)

    # bind_norms=true: all K equal
    prob2   = multi_dataset_problem(base, datasets...; bind_norms=true)
    result2 = fit(prob2, LevenbergMarquadt())
    @test result2[1].u[1] ≈ result2[2].u[1]  rtol=1e-4
    @test result2[1].u[1] ≈ result2[3].u[1]  rtol=1e-4
end
