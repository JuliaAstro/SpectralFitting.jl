"""
    XS_PhotoelectricAbsorption(ηH)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, XS_PhotoelectricAbsorption())
```

```
             XS_PhotoelectricAbsorption
     ┌────────────────────────────────────────┐
   1 │       ...''''''''''''''''''''''''''''''│
     │      .'                                │
     │     :                                  │
     │    :'                                  │
     │    :                                   │
     │   :                                    │
     │   :                                    │
     │   :                                    │
     │  :                                     │
     │  :                                     │
     │  :                                     │
     │  :                                     │
     │  :                                     │
     │ :                                      │
   0 │.:                                      │
     └────────────────────────────────────────┘
      0                                     20
                       E (keV)
```
"""
@xspecmodel :C_phabs struct XS_PhotoelectricAbsorption{T} <: AbstractSpectralModel{
        T, Multiplicative,
    }
    "Equivalent hydrogen column (units of 10²² atoms per cm⁻²)."
    ηH::T
end
XS_PhotoelectricAbsorption(; ηH = FitParam(1.0)) = XS_PhotoelectricAbsorption(ηH)

"""
    XS_WarmAbsorption(ηH, Ew)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, XS_WarmAbsorption())
```

```
                    XS_WarmAbsorption
       ┌────────────────────────────────────────┐
     1 │':      ...''':'''''''''''''''''''''''''│
       │ :    .:'                               │
       │ :   .'                                 │
       │ :  .:                                  │
       │ :  :                                   │
       │ :  :                                   │
       │ : :                                    │
       │ : :                                    │
       │ : :                                    │
       │ : :                                    │
       │ ::                                     │
       │ ::                                     │
       │  :                                     │
       │  :                                     │
   0.2 │  :                                     │
       └────────────────────────────────────────┘
        0                                     20
                         E (keV)
```
"""
@xspecmodel :C_wndabs struct XS_WarmAbsorption{T} <: AbstractSpectralModel{T, Multiplicative}
    "Equivalent hydrogen column (units of 10²² atoms per cm⁻²)."
    ηH::T
    "Window energy (keV)."
    Ew::T
end
XS_WarmAbsorption(; ηH = FitParam(1.0), Ew = FitParam(1.0)) = XS_WarmAbsorption(ηH, Ew)

"""
    XS_NeutralHydrogenAbsorption(nH)

$(FIELDS)

# Example

```julia
energy = collect(range(0.1, 20.0, 100))
invokemodel(energy, XS_NeutralHydrogenAbsorption())
```

```
            XS_NeutralHydrogenAbsorption
     ┌────────────────────────────────────────┐
   1 │       ...''''''''''''''''''''''''''''''│
     │      .'                                │
     │     :                                  │
     │    :'                                  │
     │    :                                   │
     │   :                                    │
     │   :                                    │
     │   :                                    │
     │  :'                                    │
     │  :                                     │
     │  :                                     │
     │  :                                     │
     │  :                                     │
     │ :                                      │
   0 │.:                                      │
     └────────────────────────────────────────┘
      0                                     20
                       E (keV)
```
"""
@xspecmodel :C_tbabs struct XS_NeutralHydrogenAbsorption{T} <:
    AbstractSpectralModel{T, Multiplicative}
    "Neutral hydrogen column density (units of 10²² atoms per cm⁻²)"
    nH::T
end
XS_NeutralHydrogenAbsorption(; nH = FitParam(1.0)) = XS_NeutralHydrogenAbsorption(nH)

export XS_PhotoelectricAbsorption, XS_WarmAbsorption, XS_NeutralHydrogenAbsorption
