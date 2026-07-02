export AbstractFittingAlgorithm, LevenbergMarquadt, fit, fit!

abstract type AbstractFittingAlgorithm end

struct LevenbergMarquadt{T} <: AbstractFittingAlgorithm
    λ_inc::T
    λ_dec::T
end
LevenbergMarquadt(; λ_inc = 10.0, λ_dec = 0.1) = LevenbergMarquadt(λ_inc, λ_dec)

function _lsq_fit(
        f,
        x,
        y,
        cov,
        parameters,
        alg;
        verbose = false,
        max_iter = 1000,
        kwargs...,
    )
    return LsqFit.curve_fit(
        f,
        x,
        y,
        cov,
        get_value.(parameters);
        lower = get_lowerlimit.(parameters),
        upper = get_upperlimit.(parameters),
        lambda_increase = alg.λ_inc,
        lambda_decrease = alg.λ_dec,
        show_trace = verbose,
        maxIter = max_iter,
        kwargs...,
    )
end

function configuration(prob::FittingProblem; kwargs...)
    kw, config = _unpack_config(prob; kwargs...)
    if length(kw) > 0
        throw("Unknown keyword arguments: $(kw)")
    end
    return config
end

function fit(prob::FittingProblem, args...; kwargs...)
    method_kwargs, config = _unpack_config(prob; kwargs...)
    return fit(config, args...; method_kwargs...)
end

function fit(
        config::FittingConfig,
        alg::LevenbergMarquadt;
        verbose = false,
        max_iter = 1000,
        autodiff = supports_autodiff(config) ? :forward : :finite,
        method_kwargs...,
    )
    @assert fit_statistic(config) == ChiSquared() "Least squares only for χ2 statistics."

    total_objective = reduce(vcat, c.objective for c in config.data_cache)
    total_domain = reduce(vcat, c.model_domain for c in config.data_cache)
    total_cov = reduce(vcat, c.covariance for c in config.data_cache)

    function _invoke_wrapper(_, u)
        return vcat(calculate_objective!(config, u)...)
    end

    lsq_result = _lsq_fit(
        _invoke_wrapper,
        total_domain,
        total_objective,
        total_cov,
        config.u0,
        alg;
        verbose = verbose,
        max_iter = max_iter,
        autodiff = autodiff,
        method_kwargs...,
    )
    params = LsqFit.coef(lsq_result)
    σ = try
        LsqFit.standard_errors(lsq_result)
    catch e
        @warn "No parameter uncertainty estimation due to error: $e"
        nothing
    end

    return finalize_result(config, params, lsq_result; σparams = σ)
end

function fit(
        config::FittingConfig,
        optim_alg;
        verbose = false,
        autodiff = _determine_ad_backend(config),
        method_kwargs...,
    )
    if verbose == true
        @warn "Verbose not yet supported for these fitting algorithms."
    end

    function _stat_objective(params, _)
        out = calculate_objective!(config, params)
        m = sum(1:model_count(config.prob)) do i
            dc = config.data_cache[i]
            measure(fit_statistic(config), out[i], dc.objective, dc.variance)
        end
        return m
    end

    if !(autodiff isa Optimization.SciMLBase.NoAD) && (!supports_autodiff(config))
        error("Model does not support automatic differentiation.")
    end

    u0 = get_value.(config.u0)
    # lb, ub = _determine_bounds(config, autodiff)

    # # build problem and solve
    # # TODO: iip == false fails to solve, but this works fine?
    opt_f = Optimization.OptimizationFunction{true}(_stat_objective, autodiff)
    opt_prob = Optimization.OptimizationProblem{true}(opt_f, u0, nothing)

    sol = Optimization.solve(opt_prob, optim_alg; method_kwargs...)

    # # TODO: temporary fix for type instabilities in Optimizations.jl
    return finalize_result(config, sol.u, sol)
end

function fit!(prob::FittingProblem, args...; kwargs...)
    result = fit(prob, args...; kwargs...)
    @assert length(result.config.parameter_bindings) == 1 "Can only update model when there is a single result slice"
    update_model!(only(prob.model.m), result[1])
    return result
end

_determine_bounds(_, ::Optimization.SciMLBase.NoAD) = (nothing, nothing)
_determine_bounds(config, _) = (get_lowerlimit.(config.u0), get_upperlimit.(config.u0))

function _determine_ad_backend(config)
    if supports_autodiff(config)
        return Optimization.ADTypes.AutoForwardDiff()
    else
        return Optimization.SciMLBase.NoAD()
    end
end
