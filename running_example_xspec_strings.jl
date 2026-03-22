# Running example for Issue #187 — XSPEC model string parser
# Self-contained: only tests parse_xspec_model_string (no LibXSPEC needed)

# Inline the lookup table and parser directly
const XSPEC_MODEL_NAMES = Dict{String,Symbol}(
    "powerlaw"=>"XS_PowerLaw" |> Symbol,"po"=>:XS_PowerLaw,
    "cutoffpl"=>:XS_CutOffPowerLaw,"bbody"=>:XS_BlackBody,"bb"=>:XS_BlackBody,
    "bremss"=>:XS_BremsStrahlung,"kerrdisk"=>:XS_KerrDisk,"kyrline"=>:XS_KyrLine,
    "laor"=>:XS_Laor,"diskline"=>:XS_DiskLine,"gaussian"=>:XS_Gaussian,"gaus"=>:XS_Gaussian,
    "jet"=>:XS_Jet,"optxagnf"=>:XS_Optxagnf,"diskbb"=>:XS_DiskBlackBody,
    "phabs"=>:XS_PhotoelectricAbsorption,"wndabs"=>:XS_WarmAbsorption,
    "tbabs"=>:XS_NeutralHydrogenAbsorption,"cflux"=>:XS_CalculateFlux,"kerrconv"=>:XS_Kerrconv,
)

@enum _TK::UInt8 _NAME _STAR _PLUS _LP _RP _EOF
struct _Tok k::_TK; v::String end

function _lex(src)
    toks=_Tok[]; s=strip(lowercase(src)); i=firstindex(s)
    while i<=lastindex(s)
        c=s[i]
        if isspace(c) i=nextind(s,i)
        elseif c=='*' push!(toks,_Tok(_STAR,"")); i=nextind(s,i)
        elseif c=='+' push!(toks,_Tok(_PLUS,"")); i=nextind(s,i)
        elseif c=='(' push!(toks,_Tok(_LP,"")); i=nextind(s,i)
        elseif c==')' push!(toks,_Tok(_RP,"")); i=nextind(s,i)
        elseif isletter(c)||isdigit(c)
            j=i
            while j<=lastindex(s)&&(isletter(s[j])||isdigit(s[j])||s[j]=='_') j=nextind(s,j) end
            push!(toks,_Tok(_NAME,s[i:prevind(s,j)])); i=j
        else error("Unexpected '$c'") end
    end
    push!(toks,_Tok(_EOF,"")); toks
end

mutable struct _P toks::Vector{_Tok}; pos::Int; src::String end
_peek(p::_P)=p.toks[p.pos]
function _adv!(p::_P) t=p.toks[p.pos]; p.pos+=1; t end

function _expr!(p)
    L=_term!(p)
    while _peek(p).k==_PLUS _adv!(p); R=_term!(p); L=Expr(:call,:+,L,R) end
    L
end
function _term!(p)
    L=_fac!(p)
    while _peek(p).k==_STAR _adv!(p); R=_fac!(p); L=Expr(:call,:*,L,R) end
    L
end
function _fac!(p)
    t=_peek(p)
    if t.k==_LP _adv!(p); inner=_expr!(p)
        _peek(p).k==_RP||error("Missing ) in \"$(p.src)\""); _adv!(p); return inner
    elseif t.k==_NAME
        _adv!(p); nm=t.v
        haskey(XSPEC_MODEL_NAMES,nm)&&return Expr(:call,XSPEC_MODEL_NAMES[nm])
        error("Unknown model '$nm' in \"$(p.src)\"\n  See XSPEC_MODEL_NAMES for supported names.")
    elseif t.k==_EOF error("Unexpected end in \"$(p.src)\"")
    else error("Unexpected token in \"$(p.src)\"") end
end

function parse_xspec_model_string(s::AbstractString)
    isempty(strip(s))&&error("Empty string")
    p=_P(_lex(s),1,String(s)); e=_expr!(p)
    _peek(p).k==_EOF||error("Extra tokens in \"$s\"")
    e
end

# ─── DEMO ────────────────────────────────────────────────────────────────────
println("="^60)
println("  XSPEC Model String Demo — Issue #187")
println("="^60)

println("\n📌 Demo 1: Parsing XSPEC strings → Julia Expr\n")
cases = [
    ("powerlaw",                          "Single additive model"),
    ("phabs*powerlaw",                    "AGN: phabs=PhotoelectricAbsorption"),
    ("tbabs*powerlaw",                    "tbabs=NeutralHydrogenAbsorption"),
    ("powerlaw+bbody",                    "Two additive components"),
    ("phabs*(powerlaw+bbody)",            "Absorption on BOTH"),
    ("phabs*powerlaw+bbody",              "Absorption on powerlaw ONLY"),
    ("tbabs*(diskbb+powerlaw+gaussian)",  "BH XRB: disc+PL+Fe line"),
    ("cflux*powerlaw",                    "Convolution flux"),
    ("kerrconv*diskline",                 "Kerr broadening"),
    ("phabs*(laor+powerlaw)",             "Relativistic Fe line"),
]
for (s,desc) in cases
    e=parse_xspec_model_string(s)
    println("  \"$s\"")
    println("   → $e")
    println("   [$desc]\n")
end

println("─"^60)
println("\n📌 Demo 2: Precedence (* > +) — astrophysically critical!\n")
eA=parse_xspec_model_string("phabs*powerlaw+bbody")
eB=parse_xspec_model_string("phabs*(powerlaw+bbody)")
println("  A: \"phabs*powerlaw+bbody\"  → $eA")
println("     [bbody is UNABSORBED]")
println("  B: \"phabs*(powerlaw+bbody)\" → $eB")
println("     [BOTH components absorbed]")
@assert eA!=eB "ERROR: these should be different!"
println("\n  ✅ Correctly different — precedence works!")

println("\n─"^60)
println("\n📌 Demo 3: Aliases and case-insensitivity\n")
@assert parse_xspec_model_string("po")==parse_xspec_model_string("powerlaw")
println("  ✅ po == powerlaw")
@assert parse_xspec_model_string("gaus")==parse_xspec_model_string("gaussian")
println("  ✅ gaus == gaussian")
@assert parse_xspec_model_string("PHABS*POWERLAW")==parse_xspec_model_string("phabs*powerlaw")
println("  ✅ PHABS*POWERLAW == phabs*powerlaw")
@assert parse_xspec_model_string("TbAbs*PowerLaw")==parse_xspec_model_string("tbabs*powerlaw")
println("  ✅ TbAbs*PowerLaw == tbabs*powerlaw")

println("\n─"^60)
println("\n📌 Demo 4: Error handling\n")
for (s,reason) in [("apec","not in package"),("*powerlaw","leading *"),
                   ("powerlaw*","trailing *"),("","empty")]
    try parse_xspec_model_string(s); println("  ❌ \"$s\" should error!")
    catch e; println("  ✅ \"$s\" → $(sprint(showerror,e)[1:min(60,end)])…")
    end
end

println("\n─"^60)
println("\n📌 Demo 5: XSPEC_MODEL_NAMES ($(length(XSPEC_MODEL_NAMES)) entries)\n")
for (k,v) in sort!(collect(XSPEC_MODEL_NAMES))
    println("  \"$k\" → $v")
end

println()
println("="^60)
println("  ✅  ALL DEMOS PASSED — implementation correct!")
println("="^60)
