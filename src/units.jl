module SpectralUnits

using Reexport
@reexport using Unitful

Unitful.register(SpectralUnits)

@dimension 𝐄 "events" Events
@refunit counts "counts" Counts 𝐄 false

const localpromotion = Unitful.promotion
function __init__()
    Unitful.register(SpectralUnits)
    merge!(Unitful.promotion, localpromotion)
    return
end

const _counts = typeof(u"counts")
const _counts_kev = typeof(u"counts / keV")
const _rate = typeof(u"counts / s")
const _rate_kev = typeof(u"counts / (s * keV)")
const RateOrCount = Union{<:_counts, <:_rate, <:_counts_kev, <:_rate_kev}

export counts

end # module
