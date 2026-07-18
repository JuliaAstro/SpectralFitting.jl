using Test, SpectralFitting, XSPECModels

testdir = get(
    ENV,
    "SF_TEST_SUITE_DATA",
    "/home/lilith/developer/jl/spectral-fitting-test-suite/sample-data",
)
@show testdir

has_test_dir = isdir(testdir)

if !has_test_dir
    @warn "No test data found. Skipping some tests."
end

if !has_test_dir && get(ENV, "CI", "false") != "false"
    error("Missing test dir on CI")
end

# download requisite model data
ALL_NEEDED_MODELS = [
    XS_PowerLaw,
    XS_BlackBody,
    XS_BremsStrahlung,
    XS_Laor,
    XS_DiskLine,
    XS_Gaussian,
    XS_PhotoelectricAbsorption,
    XS_WarmAbsorption,
    PhotoelectricAbsorption,
]

for M in ALL_NEEDED_MODELS
    SpectralFitting.download_model_data(M)
end

@testset "Main" begin
    include("test-general-xspec-models.jl")
    include("test-xspec.jl")
    include("test-xspec-models.jl")
end

@testset "XSPEC model strings (Issue #187)" begin
    include("test-xspec-string.jl")
end

@testset "Integration" begin
    include("test-validity.jl")
    include("test-sample-data.jl")
    include("test-fitting.jl")

    if has_test_dir
        include("test-sample-data.jl")
    else
        @warn "Skipping dataset tests."
    end
end
