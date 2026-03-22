module XSPECModels

using DocStringExtensions
using LibXSPEC_jll
using SpectralFitting
import SpectralFitting: _check_model_directory_present, register_model_data

const LIBXSPEC_STORAGE_PATH = joinpath(LibXSPEC_jll.artifact_dir, "spectral", "modelData")

include("ccall-wrapper.jl")

# include xspec models
include("additive.jl")
include("multiplicative.jl")
include("convolutional.jl")
include("xspec_string.jl")
export XSPEC_MODEL_NAMES, parse_xspec_model_string, xspec_model_string, @xspec_str

function register_xspec_data()
    push!(SpectralFitting.ALL_STORAGE_PATHS, LIBXSPEC_STORAGE_PATH)
    register_model_data(XS_KerrDisk, "kerrtable.fits"; root = LIBXSPEC_STORAGE_PATH)
    register_model_data(XS_Kerrconv, "kerrtable.fits"; root = LIBXSPEC_STORAGE_PATH)
    register_model_data(XS_KyrLine, "KBHline01.fits"; root = LIBXSPEC_STORAGE_PATH)
    register_model_data(XS_Laor, "ari.mod"; root = LIBXSPEC_STORAGE_PATH)
end


function __init__()
    # init HEASOFT
    if get(ENV, "SPECTRAL_FITTING_XSPEC_INIT", "") == ""
        # push our storage path
        register_xspec_data()
        ccall((:FNINIT, libXSFunctions), Cvoid, ())
        # set an environment variable so we don't accidentally init again
        ENV["SPECTRAL_FITTING_XSPEC_INIT"] = "true"
    end
end

end # module XSPECModels
