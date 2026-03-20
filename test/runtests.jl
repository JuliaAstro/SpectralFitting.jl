using Test, SpectralFitting

testdir = get(
    ENV,
    "SF_TEST_SUITE_DATA",
    expanduser("~/.julia/dev/spectral-fitting-test-suite/sample-data"),
)
@show testdir

has_test_dir = isdir(testdir)
if !has_test_dir
    @warn "No test data found. Skipping some tests."
end
if !has_test_dir && get(ENV, "CI", false)
    error("Missing test dir on CI")
end

SpectralFitting.download_model_data(PhotoelectricAbsorption)

include("dummies.jl")
include("utils.jl")

@testset "api" verbose = true begin
    include("models/test-model-api.jl")
end

@testset "composite-algebra" verbose = true begin
    include("composite/test-algebra.jl")
    include("composite/test-invocation.jl")
end

@testset "model-library" verbose = true begin
    include("models/test-table-models.jl")
    include("models/test-surrogate-models.jl")
    include("models/test-auto-cache.jl")
    include("models/test-as-convolution.jl")
    include("models/test-copy.jl")
    include("models/test-general-models.jl")
    include("models/test-julia-models.jl")
    include("models/test-wrappers.jl")
    include("models/test-patching.jl")
end

@testset "io" verbose = true begin
    # include("io/test-printing.jl")
    include("io/test-result-printing.jl")
    if has_test_dir
        include("datasets/test-ogip.jl")
    else
        @warn "Skipping OGIP dataset tests."
    end
    include("datasets/test-grouping.jl")
    include("datasets/test-units.jl")
    include("datasets/test-binning.jl")
    include("datasets/test-datasets.jl")
    include("io/test-remote-pathname-compression.jl")
end

@testset "fitting" verbose = true begin
    @time include("fitting/test-fit-simple-dataset.jl")
    @time include("fitting/test-cache.jl")
    @time include("fitting/test-binding.jl")
    @time include("fitting/test-results.jl")
    @time include("fitting/test-update-model.jl")

    @time include("fitting/test-fit-powerlaw.jl")
    @time include("fitting/test-models.jl")

    @time include("fitting/test-fit-multi.jl")
    @time include("fitting/test-fit-optim.jl")

    @time include("fitting/test-domain-extension.jl")
end

@testset "simulation" verbose = true begin
    include("simulation/test-simulation.jl")
    if has_test_dir
        include("simulation/test-sample-data-sim.jl")
    else
        @warn "Skipping simulating real observatory tests."
    end
end

using Aqua

# little bit of aqua
Aqua.test_undefined_exports(SpectralFitting)
Aqua.test_unbound_args(SpectralFitting)
