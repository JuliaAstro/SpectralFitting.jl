using Test
using SpectralFitting

xmm_rmf_path = joinpath(testdir, "xmm/pn.rmf")
nustar_rmf_path = joinpath(testdir, "nustar/nu60001047002A01_sr.rmf")

# test we can read in the different response matrices correctly
rmf_xmm = OGIP.read_rmf(xmm_rmf_path)
@inferred OGIP.read_rmf(xmm_rmf_path)
rmf_nustar = OGIP.read_rmf(nustar_rmf_path)

# make sure the correct number of values have been written in
@test count(!=(0), rmf_xmm.matrix) == 928061
@test count(!=(0), rmf_nustar.matrix) == 8536159

# test reading the ancillary files
xmm_arf_path = joinpath(testdir, "xmm/pn.arf")
nustar_arf_path = joinpath(testdir, "nustar/nu60001047002A01_sr.arf")

arf_xmm = OGIP.read_ancillary_response(xmm_arf_path)
@inferred OGIP.read_ancillary_response(xmm_arf_path)
arf_nustar = OGIP.read_ancillary_response(nustar_arf_path)

@test arf_xmm.bins_low ≈ rmf_xmm.bins_low atol = 1e-7
@test arf_nustar.bins_low ≈ rmf_nustar.bins_low atol = 1e-3

# test reading spectra
xmm_spec_path = joinpath(testdir, "xmm/pn_spec_grp.fits")
nustar_spec_path = joinpath(testdir, "nustar/nu60001047002A01_sr_grp_simple.pha")

spec_xmm = OGIP.read_spectrum(xmm_spec_path)
@inferred OGIP.read_spectrum(xmm_spec_path)
spec_nustar = OGIP.read_spectrum(nustar_spec_path)

@test spec_xmm.telescope_name == "XMM"
@test spec_nustar.telescope_name == "NuSTAR"

# test reading background
xmm_backgroud_path = joinpath(testdir, "xmm/pn_bg_spec.fits")
nustar_background_path = joinpath(testdir, "nustar/nu60001047002A01_bk.pha")

bg_xmm = OGIP.read_background(xmm_backgroud_path)
@test length(bg_xmm.channels) == 4096

bg_nustar = OGIP.read_background(nustar_background_path)
@inferred OGIP.read_background(nustar_background_path)
@test length(bg_xmm.channels) == 4096

# test reading the associated paths
xmm_paths = OGIP.read_paths_from_spectrum(xmm_spec_path)

@test xmm_paths[1] == xmm_backgroud_path
@test xmm_paths[2] == xmm_rmf_path
@test xmm_paths[3] == xmm_arf_path
