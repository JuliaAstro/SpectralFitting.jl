module SpectralFitting

import Base
import Printf
import Downloads
import Statistics
import Pkg.MiniProgressBars: MiniProgressBar, start_progress, end_progress, show_progress
import Distributions
import ConstructionBase

import Random
import FITSFiles

using SparseArrays
using Surrogates
using LinearAlgebra
using FileIO
using Interpolations
import DataInterpolations
using SpecialFunctions
using Polynomials
using PreallocationTools
using EnumX

# fitting backends
import ForwardDiff
import LsqFit
import Optimization

using DocStringExtensions

# non-General packages
# TODO: remove this as a dependency
using MultiLinearInterpolations

abstract type AbstractInstrument end

# unitful units
include("units.jl")
SpectralUnits.@reexport using .SpectralUnits


include("utils.jl")
include("print-utilities.jl")
include("support.jl")
include("statistics.jl")

include("fitparam.jl")
include("param-cache.jl")
include("abstract-models.jl")

include("composite-models.jl")

include("meta-models/wrappers.jl")
include("meta-models/table-models.jl")
include("meta-models/surrogate-models.jl")

include("datasets/ogip.jl")
include("datasets/datasets.jl")
include("datasets/binning.jl")
include("datasets/grouping.jl")
include("datasets/injectivedata.jl")
include("datasets/binneddata.jl")

include("model-data-io.jl")

# include fitting api
include("fitting/problem.jl")
include("fitting/config.jl")
include("fitting/result.jl")
include("fitting/methods.jl")

include("simulate.jl")
include("fitting/goodness.jl")

# include julia models
include("julia-models/model-utilities.jl")
include("julia-models/additive.jl")
include("julia-models/multiplicative.jl")
include("julia-models/convolutional.jl")

# wrapper models
include("meta-models/caching.jl")
include("meta-models/functions.jl")
include("meta-models/parameter-patch.jl")

include("plots-recipes.jl")

function __init__()
    # check if we have the minimum model data already
    _check_model_directory_present()
end

end # module
