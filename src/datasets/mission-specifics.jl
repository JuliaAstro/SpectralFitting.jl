struct NuSTAR <: AbstractInstrument end
struct XmmEPIC <: AbstractInstrument end

function NuStarData(spec_path; rmf_matrix_index = 3, rmf_energy_index = 2, kwargs...)
    return OGIPDataset(
        spec_path;
        tag = NuSTAR(),
        rmf_energy_index = rmf_energy_index,
        rmf_matrix_index = rmf_matrix_index,
        kwargs...,
    )
end

function XmmData(spec_path; rmf_matrix_index = 2, rmf_energy_index = 3, kwargs...)
    return OGIPDataset(
        spec_path;
        tag = XmmEPIC(),
        rmf_energy_index = rmf_energy_index,
        rmf_matrix_index = rmf_matrix_index,
        kwargs...,
    )
end

export NuStarData, XmmData, XmmEPIC
