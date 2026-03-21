using Printf
using RecipesBase

plotting_domain(dataset::AbstractDataset) = SpectralFitting.spectrum_energy(dataset)
plotting_domain(dataset::InjectiveData) = dataset.domain

@recipe function _plotting_func(dataset::InjectiveData; data_layout = OneToOne())
    seriestype --> :scatter
    markersize --> 1.0
    markershape --> :none
    markerstrokecolor --> :auto
    yerr -> dataset.codomain_variance
    xerr -> dataset.domain_variance
    label --> make_label(dataset)
    minorgrid --> true
    dataset.domain, dataset.codomain
end

@recipe function _plotting_func(
    dataset::AbstractDataset;
    data_layout = ContiguouslyBinned(),
)
    seriestype --> :scatter
    markersize --> 0.5
    markershape --> :none
    (rate, rateerror) = (
        make_objective(data_layout, dataset),
        make_objective_variance(data_layout, dataset),
    )
    _yerr = sqrt.(rateerror)
    yerr --> _yerr
    _xerr = SpectralFitting.bin_widths(dataset) ./ 2
    xerr --> _xerr
    markerstrokecolor --> :auto
    xlabel --> "Energy (keV)"
    ylabel --> SpectralFitting.objective_units(dataset)
    label --> SpectralFitting.make_label(dataset)
    minorgrid --> true
    x = plotting_domain(dataset)
    I = @. !isinf(x) && !isinf(rate)
    @views (x[I], rate[I])
end

# ratio plots
@userplot plotbackground
@recipe function _plotting_func(p::plotbackground)
    data = p.args[1]
    background_dataset(data)
end

@recipe _plotting_func(::Type{<:FitResult}, result::FitResult) = result[1]

@recipe function _plotting_func(slice::FitResultSlice)
    label -->
    statistic_symbol(fit_statistic(slice.parent.config)) *
    Printf.@sprintf("=%.2f, ", slice.stats) *
    reduced_statistic_symbol(fit_statistic(slice.parent.config)) *
    Printf.@sprintf("=%.2f", reduced_statistic(slice))
    seriestype --> :stepmid
    x = plotting_domain(slice)
    y = calculate_objective!(slice, slice.u)
    @assert size(x) == size(y) "$(size(x)) == $(size(y))"
    x, y
end

# ratio plots
@userplot RatioPlot
@recipe function _plotting_func(
    r::RatioPlot;
    datacolor = :auto,
    modelcolor = :auto,
    label = :auto,
)
    if length(r.args) != 1 || !(typeof(r.args[1]) <: AbstractFittingResult)
        error(
            "Ratio plots first argument must be `AbstractDataset` and second argument of type `AbstractFittingResult`.",
        )
    end

    result = r.args[1] isa FittingResult ? r.args[1][1] : r.args[1]
    data = get_dataset(result)
    x = plotting_domain(data)
    y = invoke_result(result, result.u)

    y_ratio = @. result.objective / y

    ylabel --> "Ratio [data / model]"
    xlabel --> "Energy (keV)"
    minorgrid --> true

    if (label == :auto)
        label = make_label(data)
    end

    @series begin
        linestyle --> :dash
        seriestype --> :hline
        label --> false
        color --> modelcolor
        [1.0]
    end

    @series begin
        markerstrokecolor --> datacolor
        label --> label
        seriestype --> :scatter
        markershape --> :none
        markersize --> 0.5
        yerror --> sqrt.(result.variance) ./ y
        xerror --> SpectralFitting.bin_widths(data) ./ 2
        x, y_ratio
    end
end

# residual plots
# TODO: multiple datasets require repeated calls to this function (write a wrapper later)
@userplot ResidualPlot
@recipe function _plotting_fun(r::ResidualPlot)
    slices = if r.args[1] isa FitResult
        [r.args[1][i] for i = 1:result_count(r.args[1])]
    else
        @assert r.args[1] isa FitResultSlice
        [r.args[1]]
    end

    for (i, result) in enumerate(slices)
        x = plotting_domain(get_dataset(result))
        y_residuals = residuals(result)
        @series begin
            yscale := :identity
            xlabel --> "Energy (keV)"
            link := :x
            seriestype --> :stepmid
            ylabel := "Residuals"
            label := false
            fill --> (0, 0.3, :auto)
            (x, y_residuals)
        end
    end
end

# unfolded plots
@userplot UnfoldedPlot
@recipe function _plotting_fun(
    r::UnfoldedPlot;
    pow = 2,
    datacolor = :auto,
    modelcolor = :auto,
    label = :auto,
)
    if length(r.args) != 1 || !(typeof(r.args[1]) <: AbstractFittingResult)
        error(
            "Unfolded plots first argument must be `AbstractDataset` and second argument of type `AbstractFittingResult`.",
        )
    end

    result = r.args[1] isa FittingResult ? r.args[1][1] : r.args[1]
    data = get_dataset(result)
    model = get_model(result)

    x_plot = plotting_domain(data)

    y_folded = invoke_result(result, result.u)
    y_unfolded = invokemodel(data.data, model, result.u)
    Δx = bin_widths(data)

    unfolded_spectrum = @. x_plot^pow * result.objective * (y_unfolded / y_folded) / Δx
    unfolded_std = @. sqrt(result.variance) * x_plot^pow * (y_unfolded / y_folded) / Δx
    unfolded_model = @. x_plot^pow * y_unfolded / Δx

    # TODO: use unitful units to automatically label the x and y axes - this would be a great feature
    pow_prefix = pow == 0 ? "" : "E^$pow ("
    pow_suffix = pow == 0 ? "" : ")"
    ylabel --> pow_prefix * "Photons cm^-2 s^-1 keV^-1" * pow_suffix
    xlabel --> "Energy (keV)"
    minorgrid --> true

    if (label == :auto)
        label = make_label(data)
    end

    # unfolded data
    @series begin
        markerstrokecolor --> datacolor
        label --> label
        legend --> :none
        seriestype --> :scatter
        markershape --> :none
        markersize --> 0.5
        yerror --> unfolded_std
        xerror --> SpectralFitting.bin_widths(data) ./ 2
        x_plot, unfolded_spectrum
    end

    # unfolded model
    @series begin
        markerstrokecolor --> modelcolor
        label --> label
        legend --> :none
        seriestype --> :line
        markershape --> :none
        markersize --> 0.5
        x_plot, unfolded_model
    end
end
