"""
    finite_diff_kernel!(f::Function, flux, energy)

Calculates the finite difference of the function `f` over the energy bin between the high and
low bin edges, via
```math
c_i = f(E_{i,\\text{high}}) - f(E_{i,\\text{low}}),
```
similar to evaluating the limits of the integral between ``E_{i,\\text{high}}`` and ``E_{i,\\text{low}}``.

This utility function is primarily used for [`Additive`](@ref) models to ensure the flux per
bin is normalised for the energy over the bin.
"""
@inline @fastmath function finite_diff_kernel!(f::Function, flux, energy)
    E1 = f(first(energy))
    @inbounds for i in eachindex(flux)
        E2 = f(energy[i + 1])
        flux[i] = E2 - E1
        E1 = E2
    end
end

@inline @fastmath function integration_kernel!(f::Function, flux, energy)
    @inbounds for i in eachindex(flux)
        E = energy[i]
        δE = energy[i + 1] - E
        flux[i] = f(E, δE)
    end
end
