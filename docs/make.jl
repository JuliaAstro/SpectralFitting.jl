ENV["GKSwstype"] = "100"

using Documenter
using Documenter.Remotes: GitHub
using SpectralFitting, XSPECModels

SpectralFitting.download_all_model_data()

makedocs(
    warnonly = [:cross_references, :autodocs_block, :missing_docs],
    modules = [SpectralFitting, XSPECModels],
    clean = true,
    sitename = "SpectralFitting.jl",
    authors = "Fergus Baker <fergus@cosroe.com>",
    repo = GitHub("JuliaAstro/SpectralFitting.jl"),
    pages = [
        "Home" => "index.md",
        "Walkthrough" => "walkthrough.md",
        "Examples" => [
            "Diverse examples" => "examples/examples.md",
            "A quick guide" => "examples/sherpa-example.md",
        ],
        "Fitting" => [
            "parameters.md",
            # "fitting.md",
        ],
        "Models" => [
            "Using models" => "models/using-models.md",
            "Model index" => "models/models.md",
            "Composite models" => "models/composite-models.md",
            "Wrapper models" => "models/wrapper-models.md",
            "Surrogate models" => "models/surrogate-models.md",
            "XSPEC models" => "models/xspec-models.md",
        ],
         "Transitioning from XSPEC" => "transitioning-from-xspec.md", 
        "Datasets" => [
            "Using datasets" => "datasets/datasets.md",
            "Mission support" => "datasets/mission-support.md",
        ],
        "Why & How" => "why-and-how.md",
        "Reference" => "reference.md",
    ],
)

deploydocs(;
    repo = "github.com/JuliaAstro/SpectralFitting.jl.git",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#"], # Restrict to minor releases
)
