"""
    PowerLaw

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, PowerLaw())
```

```
                        PowerLaw
       ┌────────────────────────────────────────┐
   0.5 │                                        │
       │:                                       │
       │:                                       │
       │:                                       │
       │:                                       │
       │:                                       │
       │:                                       │
       │ :                                      │
       │ :                                      │
       │  :                                     │
       │   :.                                   │
       │    ':..                                │
       │        ''':......                      │
       │                  ''''''''''''''........│
     0 │                                        │
       └────────────────────────────────────────┘
        0                                     20
                         E (keV)
```
"""
struct PowerLaw{T} <: AbstractSpectralModel{T, Additive}
    "Normalisation."
    K::T
    "Photon index."
    a::T
end
PowerLaw(; K = FitParam(1.0), a = FitParam(2.0)) = PowerLaw(K, a)
@inline @fastmath function invoke!(flux, energy, model::PowerLaw)
    let a = model.a
        α = 1 - a
        if α ≈ 0.0
            finite_diff_kernel!(log, flux, energy)
        else
            α⁻¹ = inv(α)
            finite_diff_kernel!(flux, energy) do E
                α⁻¹ * E^α
            end
        end
    end
end

"""
    BlackBody

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, BlackBody())
```

```
                        BlackBody
       ┌────────────────────────────────────────┐
   0.2 │                                        │
       │                                        │
       │                                        │
       │                                        │
       │                                        │
       │                                        │
       │      .:''':..                          │
       │     :'      '':.                       │
       │   .'           ':.                     │
       │  .:               '..                  │
       │  :                  ':.                │
       │ .'                     ':..            │
       │ :                         ''...        │
       │:                              '''....  │
     0 │:                                    '''│
       └────────────────────────────────────────┘
        0                                     20
                         E (keV)
```
"""
struct BlackBody{T} <: AbstractSpectralModel{T, Additive}
    "Normalisation."
    K::T
    "Temperature (keV)."
    kT::T
end
BlackBody(; K = FitParam(1.0), kT = FitParam(3.0)) = BlackBody{typeof(K)}(K, kT)
@inline function invoke!(flux, energy, model::BlackBody)
    return let kT = model.kT
        integration_kernel!(flux, energy) do E, δE
            8.0525 * E^2 * δE / (kT^4 * (exp(E / kT) - 1))
        end
    end
end

"""
    BremsStrahlung

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, BremsStrahlung())
```

```
                   BremsStrahlung
     ┌────────────────────────────────────────┐
   3 │.                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │:                                       │
     │ :                                      │
   0 │ ':.....................................│
     └────────────────────────────────────────┘
      0                                     20
                       E (keV)
```
"""
struct BremsStrahlung{T} <: AbstractSpectralModel{T, Additive}
    "Normalisation."
    K::T
    "Temperature (keV)."
    kT::T
    "Helium to hydrogen ratio"
    ab::T
end
function BremsStrahlung(; K = FitParam(1.0), kT = FitParam(7.0), ab = FitParam(0.085))
    return BremsStrahlung{typeof(K)}(K, kT, ab)
end
@inline function invoke!(flux, energy, model::BremsStrahlung)
    return let kT = model.kT, ab = model.ab
        integration_kernel!(flux, energy) do E, δE
            γ, born_approx = _Internal_BremsStrahlung.born_approx(E, kT, 1)

            g1 =
                _Internal_BremsStrahlung.born_approx_sufficient(γ) ? born_approx :
                _Internal_BremsStrahlung.gaunt(E, kT, γ, born_approx)
            g2 =
                _Internal_BremsStrahlung.born_approx_sufficient(4γ) ? born_approx :
                _Internal_BremsStrahlung.gaunt(E, kT, 4γ, born_approx)

            gaunt = g1 + 4 * ab * g2
            gaunt * exp(-E / kT) * δE / (E * sqrt(kT))
        end
    end
end

# scope these implementation specifics in a module so they don't collide in the
# namespace
module _Internal_BremsStrahlung

    using ..Polynomials

    const A1 = [
        1.001; 1.004; 1.017; 1.036; 1.056; 1.121;; 1.001; 1.005; 1.017; 1.046; 1.073; 1.115;; 0.9991; 1.005; 1.03; 1.055; 1.102; 1.176;; 0.997; 1.005; 1.035; 1.069; 1.134; 1.186;; 0.9962; 1.004; 1.042; 1.1; 1.193; 1.306;; 0.9874; 0.9962; 1.047; 1.156; 1.327; 1.485;; 0.9681; 0.9755; 1.020; 1.208; 1.525; 1.965;;;
        0.3029; 0.1616; 0.04757; 0.013; 0.0049; -0.0032;; 0.4905; 0.2155; 0.08357; 0.02041; 0.00739; 0.00029;; 0.654; 0.2833; 0.08057; 0.03257; 0.00759; -0.00151;; 1.029; 0.391; 0.1266; 0.05149; 0.01274; 0.00324;; 0.9569; 0.4891; 0.1764; 0.05914; 0.01407; -0.00024;; 1.236; 0.7579; 0.326; 0.1077; 0.028; 0.00548;; 1.327; 1.017; 0.6017; 0.205; 0.0605; 0.00187;;;
        -1.323; -0.254; -0.01571; -0.001; -0.000184; 0.00008;; -4.762; -0.3386; -0.03571; -0.001786; -0.0003; 0.00001;; -6.349; -0.4206; -0.02571; -0.003429; -0.000234; 0.00005;; -13.231; -0.59; -0.04571; -0.005714; -0.000445; -0.00004;; -7.672; -0.6852; -0.0643; -0.005857; -0.00042; 0.00004;; -7.143; -0.9947; -0.12; -0.01007; -0.000851; -0.00004;; -3.175; -1.116; -0.2270; -0.01821; -0.001729; 0.00023
    ]

    const γ2 = [0.7783, 1.2217, 2.6234, 4.3766, 20.0, 70.0]
    const γ3 = [1.0, 1.7783, 3.0, 5.6234, 10.0, 30.0]
    const Alo1 = Polynomial(
        (
            -0.57721566,
            0.4227842,
            0.23069756,
            0.0348859,
            0.00262698,
            0.0001075,
            0.0000074,
        )
    )
    const Alo2 =
        Polynomial((1.0, 3.5156229, 3.089942, 1.2067492, 0.2659732, 0.0360768, 0.0045813))
    const Ahi1 = Polynomial(
        (
            1.25331414,
            0.07832358,
            0.02189568,
            0.01062446,
            0.00587872,
            0.0025154,
            0.00053208,
        )
    )

    poly3(weights, x) = weights[1] + weights[2] * x^1 + weights[3] * x^2

    function born_approx(E, kT, z::Int)
        if kT == 0 || E > 50 * kT || E == 0
            return zero(promote_type(typeof(kT), typeof(E)))
        end

        γ = 0.01358 * z^2 / kT
        if γ > 0.1
            return kurucz(E / kT, γ)
        end

        #  Calculate Born approximation
        u = E / 4kT
        bfac = u <= 1 ? Alo1(u^2) - log(u) * Alo2((4u / 7.5)^2) :
                Ahi1(-1 / u) / exp(2u) / sqrt(2u)
        born = 0.5513 * exp(2u) * bfac
        return (γ, born)
    end

    born_approx_sufficient(γ) = min(1000 * γ, 100) < 1

    function gaunt(E, kT, γ, born)
        u = E / 4kT
        u, γ1 = max(u, 0.003), min(1000 * γ, 100.0)
        n, m = N(u), M(γ1)
        p = (γ1 - γ3[m]) / γ2[m]
        G1 = @views A1[n, m, 1:3]
        G2 = @views A1[n, m + 1, 1:3]
        return ((1.0 - p) * poly3(G1, u) + p * poly3(G2, u)) * born
    end

    N(x) = x <= 0.03 ? 1 : x <= 0.3 ? 2 : x <= 1.0 ? 3 : x <= 5.0 ? 4 : x <= 15.0 ? 5 : 6
    M(x) = x <= 1.773 ? 1 : x <= 3.0 ? 2 : x <= 5.6234 ? 3 : x <= 10.0 ? 4 : x <= 30.0 ? 5 : 6

    const A2 = [
        5.40; 5.25; 5.00; 4.69; 4.48; 4.16; 3.85;;
        4.77; 4.63; 4.40; 4.13; 3.87; 3.52; 3.27;;
        4.15; 4.02; 3.80; 3.57; 3.27; 2.98; 2.70;;
        3.54; 3.41; 3.22; 2.97; 2.70; 2.45; 2.20;;
        2.94; 2.81; 2.65; 2.44; 2.21; 2.01; 1.81;;
        2.41; 2.32; 2.19; 2.02; 1.84; 1.67; 1.50;;
        1.95; 1.90; 1.80; 1.68; 1.52; 1.41; 1.30;;
        1.55; 1.56; 1.51; 1.42; 1.33; 1.25; 1.17;;
        1.17; 1.30; 1.32; 1.30; 1.20; 1.15; 1.11;;
        0.86; 1.00; 1.15; 1.18; 1.15; 1.11; 1.08;;
        0.59; 0.76; 0.97; 1.09; 1.13; 1.10; 1.08;;
        0.38; 0.53; 0.76; 0.96; 1.08; 1.09; 1.09
    ]

    tkur(γ, j::Int) = 2log10(γ) + 3 - j
    ukur(μ, k::Int) = 2log10(μ) + 9 - k

    function kurucz(μ, γ)
        #  Algorithm for low kT
        j = round(Integer, tkur(γ, 0))
        k = max(1, round(Integer, ukur(μ, 0)))
        t, u = tkur(γ, j), ukur(μ, k)
        return j >= 7 || j < 1 || k >= 12 ? zero(T) :
            (1 - t) * (1 - u) * A2[j, k] +
            t * (1 - u) * A2[j + 1, k] +
            (1 - t) * u * A2[j, k + 1] +
            t * u * A2[j + 1, k + 1]
    end

end # module _Internal_BremsStrahlung

"""
    GaussianLine

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, GaussianLine())
```

```
                       GaussianLine
        ┌────────────────────────────────────────┐
   0.09 │                                        │
        │            ..                          │
        │           .':                          │
        │           : :                          │
        │           : '.                         │
        │          .'  :                         │
        │          :   :                         │
        │          :   :                         │
        │          :   '.                        │
        │         :     :                        │
        │         :     :                        │
        │         :     :                        │
        │        .'      :                       │
        │        :       '.                      │
      0 │.......:         :......................│
        └────────────────────────────────────────┘
         0                                     20
                          E (keV)
```
"""
struct GaussianLine{T} <: AbstractSpectralModel{T, Additive}
    "Normalisation."
    K::T
    "Mean"
    μ::T
    "Standard deviation"
    σ::T
end
function GaussianLine(; K = FitParam(1.0), μ = FitParam(6.4), σ = FitParam(1.0))
    return GaussianLine{typeof(K)}(K, μ, σ)
end
@inline function invoke!(flux, energy, model::GaussianLine)
    return let μ = model.μ, σ = model.σ
        integration_kernel!(flux, energy) do Elow, δE
            E = Elow + δE / 2
            δE * inv(σ * √(2π)) * exp(-1 * (E - μ)^2 / (2 * σ^2))
        end
    end
end

"""
    DeltaLine

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, DeltaLine())
```

```
                        DeltaLine
       ┌────────────────────────────────────────┐
   0.4 │         .                              │
       │         :                              │
       │         :                              │
       │         :                              │
       │         :                              │
       │         :                              │
       │         ::                             │
       │         ::                             │
       │         ::                             │
       │         ::                             │
       │         ::                             │
       │         ::                             │
       │         ::                             │
       │         ::                             │
     0 │.........::.............................│
       └────────────────────────────────────────┘
        0                                     20
                         E (keV)
```

!!! note
    The `DeltaLine` model is not a true delta function, as this would be
    extremely difficult to define in a numerical model that needs to be able to
    propagate gradients. Instead, it is a very narrow [`GaussianLine`](@ref)
    model.
"""
struct DeltaLine{W <: Number, T} <: AbstractSpectralModel{T, Additive}
    _width::W
    "Normalisation."
    K::T
    "Energy at which the delta function spikes."
    E::T
end

function DeltaLine(; K = FitParam(1.0), E = FitParam(5.0), width = 1.0e-2)
    return DeltaLine{typeof(width), typeof(K)}(width, K, E)
end

@inline function invoke!(flux, energy, model::DeltaLine{T}) where {T}
    # we can't actually have a diract delta because that would ruin
    # the ability to run dual numbers through the system. What we can do instead
    # is have a miniscule gaussian
    return invoke!(flux, energy, GaussianLine(promote(model.K, model.E, model._width)...))
end

module _Xillver_Data

    using ..SpectralFitting

    const _XILLVER_TABLE_TYPE = TableModelInterpolation{Float32, 5}
    # backing data singleton
    _xillver_data::Union{Nothing, _XILLVER_TABLE_TYPE} = nothing

    function _read_table_interpolation(
            fits_path::Union{Nothing, <:AbstractString},
        )::_XILLVER_TABLE_TYPE
        if !isnothing(fits_path)
            return TableModelInterpolation(TableModelData(Val(5), default_path; T = Float32))
        end

        default_path = first(SpectralFitting.get_model_data_paths(XillverD5))
        global _xillver_data
        if !isnothing(_xillver_data)
            return _xillver_data
        end
        _xillver_data =
            TableModelInterpolation(TableModelData(Val(5), default_path; T = Float32))
        return _xillver_data
    end

end # module

struct XillverD5{T} <: AbstractTableModel{T, Additive}
    table::_Xillver_Data._XILLVER_TABLE_TYPE
    K::T
    Γ::T
    A_Fe::T
    logXi::T
    density::T
    inclination::T
end
register_model_data(
    XillverD5,
    (
        "https://sites.srl.caltech.edu/~javier/xillver/tables/xillverD-5.fits",
        "xillverD-5.fits",
    ),
)

function XillverD5(;
        K = FitParam(1.0),
        Γ = FitParam(2.0),
        A_Fe = FitParam(1.0),
        logXi = FitParam(2.0),
        density = FitParam(17.0),
        inclination = FitParam(45.0),
        table_path = nothing,
    )
    table = _Xillver_Data._read_table_interpolation(table_path)
    return XillverD5(table, K, Γ, A_Fe, logXi, density, inclination)
end

function invoke!(output, domain, model::XillverD5)
    spectra = SpectralFitting.interpolate_table!(
        model.table,
        model.Γ,
        model.A_Fe,
        model.logXi,
        model.density,
        model.inclination,
    )
    rebin_if_different_domains!(output, domain, model.table.data.energy_bins, spectra)
    return output
end

export PowerLaw, BlackBody, BremsStrahlung, GaussianLine, DeltaLine, XillverD5
