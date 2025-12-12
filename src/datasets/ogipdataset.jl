struct OGIPData end

mutable struct OGIPMetadata{H}
    paths::SpectralDataPaths
    observation_id::String
    exposure_id::String
    object::String
    header::H
end

function ogip_dataset(
    spec_path;
    hdu = 2,
    background = nothing,
    response = nothing,
    ancillary = nothing,
    kwargs...,
)
    paths = SpectralDataPaths(
        spec_path;
        background = background,
        response = response,
        ancillary = ancillary,
    )
    header = read_fits_header(paths.spectrum; hdu = hdu)

    obs_id = get(header, "OBS_ID", "[no observation id]")
    exposure_id = get(header, "EXP_ID", "[no exposure id]")
    object = get(header, "OBJECT", "[no object]")

    (; paths, obs_id, exposure_id, object, header)
end

function OGIPDataset(spec_path; tag = OGIPData(), kwargs...)
    info = ogip_dataset(spec_path; kwargs...)
    metadata =
        OGIPMetadata(info.paths, info.obs_id, info.exposure_id, info.object, info.header)
    SpectralData(info.paths, tag = tag, user_data = metadata)
end

make_label(data::SpectralData{T,<:Union{<:AbstractInstrument,<:OGIPData}}) where {T} =
    data.user_data.observation_id

function _printinfo(
    io,
    data::SpectralData{T,K},
) where {T,K<:Union{<:AbstractInstrument,<:OGIPData}}
    descr = """SpectralDataset{$K}:
      . Object              : $(data.user_data.object)
      . Observation ID      : $(data.user_data.observation_id)
      . Exposure ID         : $(data.user_data.exposure_id)
    """
    print(io, descr)
    print_spectral_data_info(io, data)
end

export OGIPDataset
