abstract type AbstractStatistic end
statistic_symbol(s::AbstractStatistic) = Base.typename(typeof(s)).name
reduced_statistic_symbol(s::AbstractStatistic) = statistic_symbol(s) * "_reduced"

struct ChiSquared <: AbstractStatistic end

statistic_symbol(::ChiSquared) = "χ²"
reduced_statistic_symbol(::ChiSquared) = "χᵥ²"

measure(::ChiSquared, y, ŷ, σ²) = sum(@.((y - ŷ)^2 / σ²))

struct Cash <: AbstractStatistic end

statistic_symbol(::Cash) = "C-stat"

function measure(::Cash, S, ŷ, _) # ignores variance
    # avoid errors in log
    any(<(0), ŷ) && return Inf
    any(<(0), S) && return Inf
    2 * sum(@.(ŷ - S + S * (log(S) - log(ŷ))))
end

export ChiSquared, Cash, measure
