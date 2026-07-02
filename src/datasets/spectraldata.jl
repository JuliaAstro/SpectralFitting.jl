struct SpectralDataPaths
    spectrum::Union{Nothing, String}
    background::Union{Nothing, String}
    response::Union{Nothing, String}
    ancillary::Union{Nothing, String}
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(paths::SpectralDataPaths))
    descr = """SpectralFilePaths:
      . Spectrum       : $(paths.spectrum)
      . Response       : $(paths.response)
      . Background     : $(paths.background)
      . Ancillary      : $(paths.ancillary)
    """
    print(io, descr)
    return
end

function SpectralDataPaths(;
        spectrum = nothing,
        background = nothing,
        response = nothing,
        ancillary = nothing,
    )
    return SpectralDataPaths(spectrum, background, response, ancillary)
end

function SpectralDataPaths(
        spec_path;
        background = nothing,
        response = nothing,
        ancillary = nothing,
    )
    background_path, response_path, ancillary_path =
        OGIP.read_paths_from_spectrum(spec_path)
    return SpectralDataPaths(
        spectrum = spec_path,
        background = isnothing(background) ? background_path : background,
        response = isnothing(response) ? response_path : response,
        ancillary = isnothing(ancillary) ? ancillary_path : ancillary,
    )
end

"""
    SpectralData{T} <: AbstractDataset


A general spectral data structure, minimally with a [`Spectrum`](@ref) and
[`ResponseMatrix`](@ref). Optionally also includes the
[`AncillaryResponse`](@ref) and a background [`Spectrum`](@ref).

## Methods

The following methods are made available through the `SpectralData`:

- [`regroup!`](@ref)
- [`restrict_domain!`](@ref)
- [`mask_energies!`](@ref)
- [`drop_channels!`](@ref)
- [`drop_bad_channels!`](@ref)
- [`drop_negative_channels!`](@ref)
- [`normalize!`](@ref)
- [`objective_units`](@ref)
- [`spectrum_energy`](@ref)
- [`bin_widths`](@ref)
- [`subtract_background!`](@ref)
- [`set_domain!`](@ref)
- [`error_statistic`](@ref)

## Constructors

The available constructors are:

    SpectralData(paths::SpectralDataPaths; kwargs...)

Using [`SpectralDataPaths`](@ref).

If the spectrum and response matrix have already been loaded separately, use

    SpectralData(
        spectrum::Spectrum,
        response::ResponseMatrix;
        # try to match the domains of the response matrix to the data
        match_domains = true,
        background = nothing,
        ancillary = nothing,
    )
"""
mutable struct SpectralData{T, Tag, D} <: AbstractDataset
    "Observed spectrum to be fitted."
    spectrum::Spectrum{T}
    "Instrument response."
    response::ResponseMatrix{T}
    "Background is optional."
    background::Union{Nothing, Spectrum{T}}
    "Ancillary response is also optional, as it may have already been folded
    through the response matrix"
    ancillary::Union{Nothing, AncillaryResponse{T}}
    "Energy of the spectrum translated from the response matrix."
    energy_low::Vector{T}
    energy_high::Vector{T}
    "Model domain to evaluate the model on"
    domain::Vector{T}
    "Mask representing which bins are to be included in the fit."
    data_mask::BitVector
    user_data::D
end

# constructor

function SpectralData(
        paths::SpectralDataPaths;
        tag = nothing,
        user_data = nothing,
        kwargs...,
    )
    spec, resp, background, ancillary = _read_all_ogip(paths; kwargs...)
    return SpectralData(spec, resp; background, ancillary, tag, user_data)
end

function SpectralData(
        spectrum::Spectrum{T},
        response::ResponseMatrix;
        # try to match the domains of the response matrix to the data
        match_domains = true,
        background = nothing,
        ancillary = nothing,
        tag = nothing,
        user_data = nothing,
    ) where {T}
    domain = _make_domain_vector(spectrum, response)
    energy_low, energy_high = _make_energy_vector(spectrum, response)
    data_mask = BitVector(fill(true, size(spectrum.data)))
    data = SpectralData{T, typeof(tag), typeof(user_data)}(
        spectrum,
        response,
        background,
        ancillary,
        energy_low,
        energy_high,
        domain,
        data_mask,
        user_data,
    )
    if !check_domains(data)
        if match_domains
            match_domains!(data)
        else
            @warn "The spectrum and response domains are unmatched. Use `match_domains!` to remedy. Results will assume you know what you're doing"
        end
    end
    return data
end

function Base.copy(data::SpectralData)
    return typeof(data)(
        data.spectrum,
        data.response,
        data.background,
        data.ancillary,
        copy(data.energy_low),
        copy(data.energy_high),
        copy(data.domain),
        copy(data.data_mask),
        copy(data.user_data),
    )
end

supports(::Type{<:SpectralData}) = (ContiguouslyBinned(),)

function _objective_to_units(dataset::SpectralData, obj, units)
    adj = if !isnothing(units)
        alt = copy(obj)
        adjust_to_units!(dataset, dataset.spectrum, alt, units)
        alt
    else
        obj
    end
    return adj[dataset.data_mask]
end

function make_objective(layout::AbstractDataLayout, dataset::SpectralData)
    obj = make_objective(layout, dataset.spectrum)
    return _objective_to_units(dataset, obj, support_units(layout))
end

function make_objective_variance(layout::AbstractDataLayout, dataset::SpectralData)
    var = make_objective_variance(layout, dataset.spectrum)
    return _objective_to_units(dataset, var, support_units(layout))
end

make_model_domain(::ContiguouslyBinned, dataset::SpectralData) = dataset.domain
make_output_domain(::ContiguouslyBinned, dataset::SpectralData) =
    folded_energy(dataset.response)

restrict_domain!(dataset::SpectralData, low, high) =
    restrict_domain!(dataset, i -> high > i > low)
mask_energies!(dataset::SpectralData, low, high) =
    mask_energies!(dataset, i -> high > i > low)

function mask_energies!(dataset::SpectralData, condition)
    J = @. !condition(dataset.energy_low) || !condition(dataset.energy_high)
    dataset.data_mask[J] .= false
    return dataset
end

function restrict_domain!(dataset::SpectralData, condition)
    mask_energies!(dataset, condition)
    return dataset
end

function _fold_transformer(T::Type, exposure_time, layout::AbstractDataLayout, R, ΔE, E)
    cache = DiffCache(construct_objective_cache(layout, T, length(E), 1))
    units = support_units(layout)
    function _transformer!!(energy, flux)
        f = rebin_if_different_domains!(get_tmp(cache, flux), E, energy, flux)
        f = R * f
        if units === u"counts / (s * keV)"
            @. f = f / ΔE
        elseif units === u"counts"
            @. f = f * exposure_time
        end
        return f
    end
    function _transformer!!(output, energy, flux)
        f = rebin_if_different_domains!(get_tmp(cache, flux), E, energy, flux)
        mul!(output, R, f)
        if units === u"counts / (s * keV)"
            @. output = output / ΔE
        elseif units === u"counts"
            @. output = output * exposure_time
        end
        return output
    end
    return _transformer!!
end

function objective_transformer(
        layout::ContiguouslyBinned,
        dataset::SpectralData{T},
    ) where {T}
    R_folded = response_matrix(dataset)
    R = R_folded[dataset.data_mask, :]
    ΔE = bin_widths(dataset)
    model_domain = response_energy(dataset.response)
    return _fold_transformer(T, dataset.spectrum.exposure_time, layout, R, ΔE, model_domain)
end

response_matrix(dataset::SpectralData) =
    has_ancillary(dataset) ? sparse(fold_ancillary(dataset.response, dataset.ancillary)) :
    dataset.response.matrix

unmasked_bin_widths(dataset::SpectralData) = dataset.energy_high .- dataset.energy_low
bin_widths(dataset::SpectralData) = unmasked_bin_widths(dataset)[dataset.data_mask]
has_background(dataset::SpectralData) = !isnothing(dataset.background)
has_ancillary(dataset::SpectralData) = !isnothing(dataset.ancillary)

function drop_bad_channels!(dataset::SpectralData)
    indices = findall(!=(GOOD_QUALITY), dataset.spectrum.quality)
    return drop_channels!(dataset, indices)
end

function drop_negative_channels!(dataset::SpectralData)
    indices = findall(<(0), dataset.spectrum.data)
    return drop_channels!(dataset, indices)
end

function drop_channels!(dataset::SpectralData, indices)
    drop_channels!(dataset.spectrum, indices)
    drop_channels!(dataset.response, indices)
    if has_background(dataset)
        drop_channels!(dataset.background, indices)
    end
    deleteat!(dataset.data_mask, indices)
    deleteat!(dataset.energy_low, indices)
    deleteat!(dataset.energy_high, indices)
    return length(indices)
end

spectrum_energy(dataset::SpectralData) =
    ((dataset.energy_low .+ dataset.energy_high) ./ 2)[dataset.data_mask]

function regroup!(dataset::SpectralData, grouping; safety_copy = false)
    grp::typeof(grouping) = if safety_copy
        copy(grouping)
    else
        grouping
    end

    itt = GroupingIterator(grp)
    for i in itt
        dataset.energy_low[i[1]] = dataset.energy_low[i[2]]
        dataset.energy_high[i[1]] = dataset.energy_high[i[3]]
    end

    if has_background(dataset)
        regroup!(dataset.background, grp)
    end
    regroup!(dataset.response, grp)
    regroup!(dataset.spectrum, grp)

    resize!(dataset.data_mask, length(itt))
    resize!(dataset.energy_low, length(itt))
    resize!(dataset.energy_high, length(itt))
    # set everything to unmasked
    dataset.data_mask .= 1
    return dataset
end

function regroup!(dataset::SpectralData; min_counts = nothing)
    if !isnothing(min_counts)
        group_min_counts!(dataset.spectrum, min_counts)
    end
    return regroup!(dataset, dataset.spectrum.grouping)
end

function normalize!(dataset::SpectralData)
    Base.depwarn("`normalize!` is deprecated. Use `set_units!` instead.", :normalize!)
    return set_units!(dataset, u"counts / (s * keV)")
end

function subtract_background!(dataset::SpectralData)
    if !has_background(dataset)
        error("No background to subtract. Did you already subtract the background?")
    end
    subtract_background!(dataset.spectrum, dataset.background)
    dataset.background = nothing
    return dataset
end

function set_domain!(dataset::SpectralData, domain)
    return dataset.domain = domain
end

objective_units(data::SpectralData) = data.spectrum.units

error_statistic(data::SpectralData) = error_statistic(data.spectrum)

preferred_units(::Type{<:SpectralData}, ::AbstractStatistic) = nothing
preferred_units(::Type{<:SpectralData}, ::ChiSquared) = u"counts / (s * keV)"
preferred_units(::Type{<:SpectralData}, ::Cash) = u"counts"

_as_unit(::Unitful.FreeUnits{U}) where {U} = first(U)

function _adjust_by_unit_difference!(
        ΔE,
        exposure_time,
        x,
        ::Unitful.FreeUnits{Units},
    ) where {Units}
    foreach(Units) do u
        if u === _as_unit(u"s")
            @. x = x * exposure_time
        elseif u === _as_unit(u"s^-1")
            @. x = x / exposure_time
        elseif u === _as_unit(u"keV")
            @. x = x * ΔE
        elseif u === _as_unit(u"keV^-1")
            @. x = x / ΔE
        else
            throw(
                ErrorException(
                    "Cannot adjust units: unhandled unit difference: $(u) (if you are expecting to be able to use specific units, please open a bug report).",
                ),
            )
        end
    end
    return nothing
end

function adjust_to_units!(data::SpectralData, s::Spectrum, x, units)
    ΔE = unmasked_bin_widths(data)
    exposure_time = s.exposure_time
    _adjust_by_unit_difference!(ΔE, exposure_time, x, units / s.units)
    return x
end

function set_units!(s::SpectralData, units)
    adjust_to_units!(s, s.spectrum, s.spectrum.data, units)
    adjust_to_units!(s, s.spectrum, s.spectrum.errors, units)
    s.spectrum.units = units
    if has_background(s)
        adjust_to_units!(s, s.background, s.background.data, units)
        adjust_to_units!(s, s.background, s.background.errors, units)
        s.background.units = units
    end
    return s
end

"""
    unfold(data::SpectralData)

Unfold the dataset using the response and ancillary, as in e.g. ISIS.

Returns a new [`Spectrum`](@abs) that has had the domain masking applied.
"""
function unfold(data::SpectralData)
    unfolded, variance = if has_ancillary(data)
        _unf = unfold(data.response, data.ancillary, data.spectrum.data)
        _var =
            !isnothing(data.spectrum.errors) ?
            unfold(data.response, data.ancillary, data.spectrum.errors) : nothing
        _unf, _var
    else
        _unf = unfold(data.response, data.spectrum.data)
        _var =
            !isnothing(data.spectrum.errors) ?
            unfold(data.response, data.spectrum.errors) : nothing
        _unf, _var
    end

    unf = remake_spectrum(data.spectrum; data = unfolded, errors = variance)
    unf.units = unf.units / u"cm^2"
    mask!(unf, data.data_mask)
    return unf
end

# internal methods

function rebin_if_different_domains!(output, data_domain, model_domain, input)
    if length(data_domain) == length(model_domain)
        @. output = input
    else
        @warn "Model and data domains differ"
        interpolated_rebin!(output, data_domain, input, model_domain)
    end
    return output
end

function _read_all_ogip(paths::SpectralDataPaths; forgiving = false)
    spec = if !isnothing(paths.spectrum)
        OGIP.read_spectrum(paths.spectrum)
    else
        if !forgiving
            @warn "No spectrum file found."
        end
        nothing
    end
    back = if !isnothing(paths.background)
        OGIP.read_background(paths.background)
    else
        if !forgiving
            @warn "No background file found."
        end
        nothing
    end
    resp = if !isnothing(paths.response)
        OGIP.read_rmf(paths.response)
    else
        if !forgiving
            throw(
                "No response file found in the header. Response must be specified with the keyword `response=PATH`.",
            )
        else
            nothing
        end
    end
    ancillary = if !isnothing(paths.ancillary)
        OGIP.read_ancillary_response(paths.ancillary)
    else
        if !forgiving
            @warn "No ancillary file found."
        end
        nothing
    end

    return (spec, resp, back, ancillary)
end

function _make_domain_vector(::Spectrum, resp::ResponseMatrix{T}) where {T}
    domain = zeros(T, length(resp.bins_low) + 1)
    # todo: check these are indeed contiguous
    domain[2:end] .= resp.bins_high
    domain[1] = resp.bins_low[1]
    return domain
end

function _make_energy_vector(spec::Spectrum, resp::ResponseMatrix{T}) where {T}
    full_domain = augmented_energy_channels(
        spec.channels,
        resp.channels,
        resp.channel_bins_low,
        resp.channel_bins_high,
    )
    high = full_domain[2:end]
    # full domain becomes the low
    resize!(full_domain, length(high))
    return full_domain, high
end

function check_domains(data::SpectralData)
    return (length(data.spectrum.channels) == length(data.response.channels)) &&
        (all(i -> i in data.response.channels, data.spectrum.channels))
end

function match_domains!(data::SpectralData)
    # drop parts of the response matrix that aren't in the spectrum
    I = filter(i -> i ∈ data.spectrum.channels, data.response.channels)
    return data.response = ResponseMatrix(
        data.response.matrix[I, :],
        data.response.channels[I],
        data.response.channel_bins_low[I],
        data.response.channel_bins_high[I],
        data.response.bins_low,
        data.response.bins_high,
    )
end

function background_dataset(data::SpectralData)
    new_data = copy(data)
    new_data.spectrum = new_data.background
    new_data.background = nothing
    return new_data
end

function rescale!(data::SpectralData)
    @. data.spectrum.data = data.spectrum.data / data.spectrum.area_scale
    @. data.spectrum.errors = data.spectrum.errors / data.spectrum.area_scale
    data.spectrum.area_scale = 1
    if has_background(data)
        rescale_background!(data)
    end
    return data
end

function rescale_background!(data::SpectralData)
    if has_background(data)
        data.background.data = _scaled_background(
            data.background.data,
            data.background.area_scale,
            data.spectrum.background_scale,
            data.background.background_scale,
        )
        data.background.errors = _scaled_background(
            data.background.errors,
            data.background.area_scale,
            data.spectrum.background_scale,
            data.background.background_scale,
        )

        data.spectrum.background_scale = 1
        data.background.background_scale = 1
        data.background.area_scale = 1
    else
        throw("No background to subtract")
    end
    return data
end

# for invoking a model on the same domain as the spectral data
function invokemodel(
        data::SpectralData,
        model::AbstractSpectralModel,
        new_free_params = nothing,
    )
    domain = make_output_domain(common_support(data, model), data)
    cache = make_parameter_cache(model)
    if !isnothing(new_free_params)
        update_free_parameters!(cache, new_free_params)
    end
    output = allocate_model_output(model, domain)
    return invokemodel!(output, domain, model, cache)[data.data_mask]
end

# printing utilities

Base.show(io::IO, @nospecialize(data::SpectralData)) =
    print(io, "SpectralData[$(data.spectrum.telescope_name)]")

_printinfo(io, data::SpectralData) = print_spectral_data_info(io, data)

function print_spectral_data_info(io, data::SpectralData{T}) where {T}
    domain = @views data.domain
    ce_min = @views prettyfloat.(minimum(data.energy_low[data.data_mask]))
    ce_max = @views prettyfloat.(maximum(data.energy_high[data.data_mask]))
    dom_min, dom_max = @views prettyfloat.(extrema(domain))

    printstyled(io, "SpectralData", color = :cyan)
    print(io, " with ")
    @views printstyled(io, length(data.energy_low[data.data_mask]), color = :cyan)
    println(io, " active channels: ")

    descr = """  . Chn. E (min/max)    : ($ce_min, $ce_max)
      . Masked channels     : $(count(==(false), data.data_mask)) / $(length(data.data_mask))
      . Model domain size   : $(length(domain))
      . Domain (min/max)    : ($dom_min, $dom_max)
    """
    print(io, descr)

    printstyled(io, "Primary Spectrum:", color = :cyan)
    println(io, "")
    _printinfo(io, data.spectrum)

    printstyled(io, "Response:", color = :cyan)
    println(io, "")
    _printinfo(io, data.response)

    if has_background(data)
        printstyled(io, "Background:", color = :cyan)
        println(io, "")
        _printinfo(io, data.background)
    else
        printstyled(io, "Background: nothing", color = :gray)
        println(io, "")
    end

    if has_ancillary(data)
        printstyled(io, "Ancillary:", color = :cyan)
        println(io, "")
        _printinfo(io, data.ancillary)
    else
        printstyled(io, "Ancillary: nothing", color = :gray)
        println(io, "")
    end
    return
end

export SpectralData,
    restrict_domain!,
    mask_energies!,
    drop_bad_channels!,
    drop_negative_channels!,
    drop_channels!,
    normalize!,
    subtract_background!,
    set_domain!,
    set_units!,
    background_dataset,
    rescale!,
    rescale_background!,
    unfold
