function goodness(
        result::FitResultSlice,
        u,
        σ;
        N = 1000,
        stat = ChiSquared(),
        distribution = Distributions.Normal,
        seed = abs(randint()),
        kwargs...,
    )
    x = similar(u)
    measures = zeros(eltype(u), N)
    config = deepcopy(result.parent.config)

    rng = Random.default_rng()
    Random.seed!(rng, seed)

    p_vector = filter!(isfree, parameter_vector(get_model(result)))

    lower_bounds = get_lowerlimit.(p_vector)
    upper_bounds = get_upperlimit.(p_vector)

    for i in eachindex(measures)
        # sample the next parameters
        for j in eachindex(u)
            mean = u[j]
            std = σ[j]
            # TODO: respect the upper and lower bounds of the parameters
            distr = Distributions.Truncated(
                distribution(mean, std),
                lower_bounds[j],
                upper_bounds[j],
            )
            x[j] = rand(distr)
        end

        simulate!(config, x; rng = rng, kwargs...)

        new_result = fit(config, LevenbergMarquadt())
        measures[i] = measure(stat, new_result)
    end

    perc = 100 * count(<(result.stats), measures) / N
    return perc, measures
end

goodness(result::FitResult, args...; kwargs...) = goodness(result[1], args...; kwargs...)

function goodness(result::FitResultSlice, σu = result.err; kwargs...)
    @assert !isnothing(σu) "σ cannot be nothing, else algorithm has no parameter intervals to sample from."
    return goodness(result, result.u, σu; kwargs...)
end

export goodness
