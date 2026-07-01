"""
    XS_PowerLaw(K, a)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, XS_PowerLaw())
```

```
                      XS_PowerLaw
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
@xspecmodel :C_powerlaw struct XS_PowerLaw{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Photon index."
    a::T
end
function XS_PowerLaw(; K = FitParam(1.0), a = FitParam(1.0))
    XS_PowerLaw(K, a)
end

"""
    XS_CutOffPowerLaw(K, Γ, Ecut, z)

$(FIELDS)

# Example

```julia
using SpectralFitting
using UnicodePlots
energy = 10 .^collect(range(-1.0, 2.0, 100))
m = invokemodel(energy, XS_CutOffPowerLaw())
lineplot(energy[1:end-1],m,xscale=:log10,yscale=:log10,xlim=(1e-1,1e2),ylim=(1e-6,1e0),xlabel="Energy (keV)",ylabel="Flux",title="XS_CutOffPowerLaw",canvas=DotCanvas)
```

```
                    XS_CutOffPowerLaw
        ┌────────────────────────────────────────┐
10⁰     │:..                                     │
        │   '':...                               │
        │        '':..                           │
        │             '''..                      │
        │                  '':..                 │
        │                      '':.              │
        │                          ':.           │
Flux    │                             '..        │
        │                               ':.      │
        │                                 '.     │
        │                                   :    │
        │                                    :.  │
        │                                     :  │
        │                                      : │
10⁻⁶    │                                       :│
        └────────────────────────────────────────┘
         10⁻¹                                 10²
                        Energy (keV)
```
"""
@xspecmodel :C_zcutoffpl struct XS_CutOffPowerLaw{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Photon index."
    Γ::T
    "Cut-off energy (keV)."
    Ecut::T
    "Redshift."
    z::T
end
function XS_CutOffPowerLaw(;
    K = FitParam(1.0),
    Γ = FitParam(2.0),
    Ecut = FitParam(15.0),
    z = FitParam(0.0, frozen = true),
)
    XS_CutOffPowerLaw(K, Γ, Ecut, z)
end

"""
    XS_BlackBody(K, T)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, XS_BlackBody())
```

```
                      XS_BlackBody
       ┌────────────────────────────────────────┐
   0.2 │                                        │
       │                                        │
       │                                        │
       │                                        │
       │                                        │
       │                                        │
       │      .:''':..                          │
       │    .:       ''.                        │
       │   .'           ':.                     │
       │   :              ''..                  │
       │  :                  ':.                │
       │ :                     '':.             │
       │.:                         ''..         │
       │:                              '':....  │
     0 │'                                    '''│
       └────────────────────────────────────────┘
        0                                     20
                         E (keV)
```
"""
@xspecmodel :C_bbody struct XS_BlackBody{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Temperature (keV)."
    T::T
end
function XS_BlackBody(; K = FitParam(1.0), T = FitParam(3.0))
    XS_BlackBody(K, T)
end

"""
    XS_BremsStrahlung(K, T)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, XS_BremsStrahlung())
```

```
                  XS_BremsStrahlung
     ┌────────────────────────────────────────┐
   2 │                                        │
     │.                                       │
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
     │'.                                      │
     │ :                                      │
   0 │  ':....................................│
     └────────────────────────────────────────┘
      0                                     20
                       E (keV)
```
"""
@xspecmodel :C_bremss struct XS_BremsStrahlung{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Plasma temperature (keV)."
    T::T
end
function XS_BremsStrahlung(; K = FitParam(1.0), T = FitParam(7.0))
    XS_BremsStrahlung(K, T)
end


"""
    XS_KerrDisk(K, lineE, index1, index2, break_r, a, θ, inner_r, outer_r)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, XS_KerrDisk())
```

```
                        XS_KerrDisk
        ┌────────────────────────────────────────┐
   0.05 │                                        │
        │                                        │
        │                                .       │
        │                               .:       │
        │                              :::       │
        │                            .:' '.      │
        │                           .:    :      │
        │                         ..'     :      │
        │                         :'      :      │
        │                       .'        :      │
        │                     .:'         :      │
        │                  .:''           :      │
        │               .::'              :      │
        │            ..:'                 :      │
      0 │.........:'''                    :......│
        └────────────────────────────────────────┘
         0                                      8
                          E (keV)
```
"""
@xspecmodel :C_kerrdisk struct XS_KerrDisk{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Rest frame line energy (keV)."
    lineE::T
    "Emissivity index for inner disk."
    index1::T
    "Emissivity index for outer disk."
    index2::T
    "Break radius separating inner and outer disk (gᵣ)."
    break_r::T
    "Dimensionless black hole spin."
    a::T
    "Disk inclination angle to line of sight (degrees)."
    θ::T
    "Inner radius of the disk in units of rₘₛ."
    inner_r::T
    "Outer radius of the disk in units of rₘₛ."
    outer_r::T
    "Redshift."
    z::T
end
function XS_KerrDisk(;
    K = FitParam(1.0),
    lineE = FitParam(6.4, frozen = true),
    index1 = FitParam(3.0, frozen = true),
    index2 = FitParam(3.0, frozen = true),
    break_r = FitParam(6.0, frozen = true),
    a = FitParam(0.998, lower_limit = 0, upper_limit = 0.998),
    θ = FitParam(30.0, lower_limit = 0, upper_limit = 90.0),
    inner_r = FitParam(1.0, frozen = true),
    outer_r = FitParam(400.0, frozen = true),
    z = FitParam(0.0, frozen = true),
)
    XS_KerrDisk(K, lineE, index1, index2, break_r, a, θ, inner_r, outer_r, z)
end


"""
    XS_KyrLine(K, a, θ_obs, inner_r, ms_flag, outer_r, lineE, α, β, break_r, z, limb)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, XS_KyrLine())
```

```
                        XS_KyrLine
        ┌────────────────────────────────────────┐
   0.05 │                                        │
        │                                        │
        │                                :       │
        │                                :.      │
        │                              :.':      │
        │                             :'  :      │
        │                            :    :      │
        │                          .'     :      │
        │                        .:'      :      │
        │                       .'        :      │
        │                     .:          :      │
        │                   .'            :      │
        │                .:'              :      │
        │            ..:'                 :      │
      0 │.........:'''                    :......│
        └────────────────────────────────────────┘
         0                                      8
                          E (keV)
```
"""
@xspecmodel :C_kyrline struct XS_KyrLine{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Dimensionless black hole spin."
    a::T
    "Observer inclination (0 is on pole, degrees)."
    θ::T
    "Inner radius of the disk in units of GM/c²"
    inner_r::T
    "0: integrate from rᵢₙ. 1: integrate from rₘₛ."
    ms_flag::T
    "Outer radius of the disk in units of GM/c²"
    outer_r::T
    "Rest frame line energy (keV)."
    lineE::T
    α::T
    β::T
    "Break radius separating inner and outer disk (GM/c²)."
    break_r::T
    "Overall Doppler shift."
    z::T
    "0: isotropic emission, 1: Laor's limb darkening, 2: Haard's limb brightening."
    limb::T
end
function XS_KyrLine(;
    K = FitParam(1.0),
    a = FitParam(0.998, lower_limit = 0, upper_limit = 1.0),
    θ = FitParam(30.0),
    inner_r = FitParam(1.0, frozen = true),
    ms_flag = FitParam(0.0, frozen = true),
    outer_r = FitParam(400.0, frozen = true),
    lineE = FitParam(6.4, frozen = true),
    α = FitParam(3.0, frozen = true),
    β = FitParam(3.0, frozen = true),
    break_r = FitParam(6.0, frozen = true),
    z = FitParam(0.0, frozen = true),
    limb = FitParam(1.0, frozen = true),
)
    XS_KyrLine(K, a, θ, inner_r, ms_flag, outer_r, lineE, α, β, break_r, z, limb)
end


"""
    XS_Laor(K, lineE, a, inner_r, outer_r, incl)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 10.0, 100))
invokemodel(energy, XS_Laor())
```

```
                          XS_Laor
        ┌────────────────────────────────────────┐
   0.06 │                                        │
        │                                        │
        │                         ::             │
        │                         ::             │
        │                        : :             │
        │                       :  :             │
        │                      :   :             │
        │                     :'   :             │
        │                   .:     :             │
        │                  :'      :             │
        │                .'        :             │
        │              .:'         :             │
        │            ..'           :             │
        │          .:'              :            │
      0 │.......:''                 :............│
        └────────────────────────────────────────┘
         0                                     10
                          E (keV)
```
"""
@xspecmodel :C_laor struct XS_Laor{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Rest frame line energy (keV)."
    lineE::T
    "Power law dependence of emissivity. Scales R⁻ᵅ."
    a::T
    "Inner radius of the accretion disk (GM/c)."
    inner_r::T
    "Outer radius of the accretion disk (GM/c)."
    outer_r::T
    "Disk inclination angle to line of sight (degrees, 0 is pole on)."
    θ::T
end
function XS_Laor(;
    K = FitParam(1.0),
    lineE = FitParam(6.4),
    a = FitParam(3.0, frozen = true),
    inner_r = FitParam(1.235, frozen = true),
    outer_r = FitParam(400.0, frozen = true),
    θ = FitParam(30.0, upper_limit = 180, frozen = true),
)
    XS_Laor(K, lineE, a, inner_r, outer_r, θ)
end

"""
    XS_DiskLine(K, lineE, β, inner_r, outer_r, incl)

$(FIELDS)

# Example

```julia
energy = collect(range(4.0, 8.0, 100))
invokemodel(energy, XS_DiskLine())
```

```
                        XS_DiskLine
        ┌────────────────────────────────────────┐
   0.09 │                                        │
        │                           .            │
        │                           :            │
        │                           ::           │
        │                         . ::           │
        │                         : ::           │
        │                         :'':           │
        │                        .'  :           │
        │                        :    :          │
        │                        :    :          │
        │                       .'    :          │
        │                       :     :          │
        │                     .:      '.         │
        │                   .:'        :         │
      0 │...............:'''           :.........│
        └────────────────────────────────────────┘
         4                                      8
                          E (keV)
```
"""
@xspecmodel :C_diskline struct XS_DiskLine{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Rest frame line energy (keV)."
    lineE::T
    "Power law dependence of emissivity. If < 10, scales Rᵅ."
    β::T
    "Inner radius of the accretion disk (GM/c)."
    inner_r::T
    "Outer radius of the accretion disk (GM/c)."
    outer_r::T
    "Disk inclination angle to line of sight (degrees, 0 is pole on)."
    θ::T
end
function XS_DiskLine(;
    K = FitParam(1.0),
    lineE = FitParam(6.7),
    β = FitParam(-2.0, frozen = true),
    inner_r = FitParam(10.0, frozen = true),
    outer_r = FitParam(1000.0, frozen = true),
    θ = FitParam(30.0, upper_limit = 180, frozen = true),
)
    XS_DiskLine(K, lineE, β, inner_r, outer_r, θ)
end

"""
    XS_Gaussian(K, E, σ)

$(FIELDS)

# Example

```julia
energy = collect(range(4.0, 8.0, 100))
invokemodel(energy, XS_Gaussian())
```

```
                        XS_Gaussian
        ┌────────────────────────────────────────┐
   0.09 │                                        │
        │            .                           │
        │           : :                          │
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
        │        :       :                       │
      0 │.......:         :......................│
        └────────────────────────────────────────┘
         0                                     20
                          E (keV)
```
"""
@xspecmodel :C_gaussian struct XS_Gaussian{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation"
    K::T
    "Line energy in keV."
    E::T
    "Line width in keV."
    σ::T
end
function XS_Gaussian(;
    K = FitParam(1.0),
    E = FitParam(6.4, frozen = true),
    σ = FitParam(1.0),
)
    XS_Gaussian(K, E, σ)
end

"""
    XS_Jet(K, mass, Dco, log_mdot, thetaobs, BulkG, phi, zdiss, B, logPrel, gmin_inj, gbreak, gmax, s1, s2, z)

$(FIELDS)

# Example

```julia
using UnicodePlots
energy = 10 .^collect(range(-8.0, 8.0, 100))
m = invokemodel(energy, XS_Jet())
lineplot(energy[1:end-1],m,xscale=:log10,yscale=:log10,xlim=(1e-8,1e8),ylim=(1e-8,1e8),xlabel="Energy (keV)",ylabel="Flux",title="XS_Jet",canvas=DotCanvas)
```

```
                        XS_Jet
        ┌────────────────────────────────────────┐
10⁸     │                                        │
        │                                        │
        │                                        │
        │                                        │
        │  :':..                                 │
        │ .'   ''.                               │
        │:        ':.                            │
Flux    │'          ':                           │
        │             '.                         │
        │              :                         │
        │               ''.....                  │
        │                     '''''..            │
        │                            ':..        │
        │                               ':.      │
10⁻⁸    │                                 ''.    │
        └────────────────────────────────────────┘
        10⁻⁸                                 10⁸
                        Energy (keV)
```
"""
@xspecmodel :C_jet struct XS_Jet{T} <: AbstractSpectralModel{T,Additive}
    "MUST BE FIXED AT UNITY as the jet spectrum normalisation is set by the relativisitic particle power."
    K::T
    "Black hole mass in solar masses."
    mass::T
    "Comoving (proper) distance in Mpc."
    Dco::T
    "log(L/L_Edd."
    log_mdot::T
    "Inclination angle (deg)."
    thetaobs::T
    "Bulk lorentz factor of the jet."
    BulkG::T
    "Angular size scale (radians) of the jet acceleration region as seen from the black hole."
    phi::T
    "Vertical distance from the black hole of the jet dissipation region (r_g)."
    zdiss::T
    "Magnetic field in the jet (Gauss)."
    B::T
    "Log of the power injected in relativisitic particles (ergs/s)."
    logPrel::T
    "Minimum lorentz factor of the injected electrons."
    gmin_inj::T
    "Lorentz factor of the break in injected electron distribution."
    gbreak::T
    "Maximum lorentz factor."
    gmax::T
    "Injected index of the electron distribution below the break."
    s1::T
    "Injected index of the electron distribution above the break."
    s2::T
    "Cosmological redshift corresponding to the comoving distance used above."
    z::T
end
function XS_Jet(;
    K = FitParam(1.0, frozen = true, lower_limit = 1.0, upper_limit = 1.0, error = 0.01),
    mass = FitParam(
        1.0e9,
        frozen = true,
        lower_limit = 1.0,
        upper_limit = 1.0e10,
        error = 1.0e7,
    ),
    Dco = FitParam(
        3350.6,
        frozen = true,
        lower_limit = 1.0,
        upper_limit = 1.0e8,
        error = 33.506,
    ),
    log_mdot = FitParam(
        -1.0,
        frozen = false,
        lower_limit = -5.0,
        upper_limit = 2.0,
        error = 0.01,
    ),
    thetaobs = FitParam(
        3.0,
        frozen = true,
        lower_limit = 0.0,
        upper_limit = 90.0,
        error = 0.03,
    ),
    BulkG = FitParam(
        13.0,
        frozen = true,
        lower_limit = 1.0,
        upper_limit = 100.0,
        error = 0.13,
    ),
    phi = FitParam(
        0.1,
        frozen = true,
        lower_limit = 0.01,
        upper_limit = 100.0,
        error = 0.001,
    ),
    zdiss = FitParam(
        1275.0,
        frozen = true,
        lower_limit = 10.0,
        upper_limit = 10000.0,
        error = 0.026,
    ),
    B = FitParam(2.6, frozen = true, lower_limit = 0.01, upper_limit = 15.0, error = 0.026),
    logPrel = FitParam(
        43.3,
        frozen = true,
        lower_limit = 40.0,
        upper_limit = 48.0,
        error = 0.433,
    ),
    gmin_inj = FitParam(
        1.0,
        frozen = true,
        lower_limit = 1.0,
        upper_limit = 1000.0,
        error = 0.01,
    ),
    gbreak = FitParam(
        300.0,
        frozen = true,
        lower_limit = 10.0,
        upper_limit = 10_000.0,
        error = 3.0,
    ),
    gmax = FitParam(
        3000.0,
        frozen = true,
        lower_limit = 1000.0,
        upper_limit = 1.0e6,
        error = 30.0,
    ),
    s1 = FitParam(1.0, frozen = true, lower_limit = -1.0, upper_limit = 1.0, error = 0.01),
    s2 = FitParam(2.7, frozen = true, lower_limit = 1.0, upper_limit = 5.0, error = 0.027),
    z = FitParam(0.0, frozen = true, lower_limit = 0.0, upper_limit = 10.0, error = 1.0),
)
    XS_Jet(
        K,
        mass,
        Dco,
        log_mdot,
        thetaobs,
        BulkG,
        phi,
        zdiss,
        B,
        logPrel,
        gmin_inj,
        gbreak,
        gmax,
        s1,
        s2,
        z,
    )
end

"""
    XS_Optxagnf()

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, XS_Optxagnf())
```
```
                                        XS_Optxagnf
                        ┌────────────────────────────────────────┐
                      0 │                                        │
                        │:                                       │
                        │:                                       │
                        │:                                       │
                        │'.                                      │
                        │ :                                      │
                        │ '.                                     │
   Flux (log scale)     │  :.                                    │
                        │   :                                    │
                        │    '.                                  │
                        │     ''':.....                          │
                        │              '''''''':.............    │
                        │                                    ''''│
                        │                                        │
                    -30 │                                        │
                        └────────────────────────────────────────┘
                         0                                     20
                                          E (keV)
```
"""
@xspecmodel :C_optxagnf struct XS_Optxagnf{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation must be frozen."
    K::T
    "Black hole mass in solar masses."
    mass::T
    "Comoving (proper) distance in Mpc."
    Dco::T
    "The Eddington ratio. log(L/L_Edd)."
    logLoLEdd::T
    "The dimensionless black hole spin."
    astar::T
    "The coronal radius in Rg = GM/c^2."
    rcor::T
    "The log of the outer radius of the disk in units of Rg."
    logrout::T
    "The electron temperature for the soft Comptonisation component (soft excess), in keV."
    kT_e::T
    "The optical depth of the soft Comptonisation component."
    τ::T
    "The spectral index of the hard Comptonisation component (‘power law’) which has temperature fixed to 100 keV."
    Gamma::T
    "the fraction of the power below rcor which is emitted in the hard comptonisation component."
    fpl::T
    "Redshift."
    z::T
end
function XS_Optxagnf(;
    K = FitParam(1.0, frozen = true, lower_limit = 1.0, upper_limit = 1.0, error = 0.01),
    mass = FitParam(
        1.0e7,
        frozen = true,
        lower_limit = 1.0,
        upper_limit = 1.0e10,
        error = 1.0e5,
    ),
    Dco = FitParam(
        100.0,
        frozen = true,
        lower_limit = 0.01,
        upper_limit = 1.0e9,
        error = 10,
    ),
    logLoLEdd = FitParam(
        -1.0,
        frozen = false,
        lower_limit = -10.0,
        upper_limit = 2.0,
        error = 0.01,
    ),
    astar = FitParam(
        0.0,
        frozen = true,
        lower_limit = 0.0,
        upper_limit = 0.998,
        error = 0.01,
    ),
    rcor = FitParam(
        10.0,
        frozen = false,
        lower_limit = 1.0,
        upper_limit = 100.0,
        error = 0.1,
    ),
    logrout = FitParam(
        5.0,
        frozen = true,
        lower_limit = 3.0,
        upper_limit = 7.0,
        error = 0.01,
    ),
    kT_e = FitParam(
        0.2,
        frozen = false,
        lower_limit = 0.01,
        upper_limit = 10.0,
        error = 0.01,
    ),
    τ = FitParam(10.0, frozen = false, lower_limit = 0.1, upper_limit = 100.0, error = 0.1),
    Gamma = FitParam(
        2.1,
        frozen = false,
        lower_limit = 1.05,
        upper_limit = 5.0,
        error = 0.01,
    ),
    fpl = FitParam(
        1.0e-4,
        frozen = false,
        lower_limit = 0.0,
        upper_limit = 1.0,
        error = 0.01,
    ),
    z = FitParam(0.0, frozen = true, lower_limit = 0.0, upper_limit = 10.0, error = 1.0),
)
    XS_Optxagnf(K, mass, Dco, logLoLEdd, astar, rcor, logrout, kT_e, τ, Gamma, fpl, z)
end

"""
- https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelDiskbb.html
"""
@xspecmodel :C_diskbb struct XS_DiskBlackBody{T} <: AbstractSpectralModel{T,Additive}
    "Normalisation."
    K::T
    "Inner disc temperature in keV."
    T::T
end
function XS_DiskBlackBody(; K = FitParam(1.0), T = FitParam(0.25))
    XS_DiskBlackBody(K, T)
end

export XS_PowerLaw,
    XS_CutOffPowerLaw,
    XS_BlackBody,
    XS_BremsStrahlung,
    XS_Laor,
    XS_DiskLine,
    XS_KerrDisk,
    XS_KyrLine,
    XS_Gaussian,
    XS_Jet,
    XS_Optxagnf,
    XS_DiskBlackBody
