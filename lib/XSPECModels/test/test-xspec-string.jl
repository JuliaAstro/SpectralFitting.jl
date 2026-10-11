# =============================================================================
# test/test-xspec-string.jl
#
# Tests for XSPEC model string parser/serialiser — Issue #187
#
# Style matches the existing XSPECModels test files.
# These tests cover ONLY the parser logic; they do NOT invoke actual XSPEC
# models (no LibXSPEC_jll needed), so they run on CI without binary deps.
# =============================================================================

using Test
using SpectralFitting
using XSPECModels

# Helper to build a :call Expr node (mirrors what parse_xspec_model_string returns)
_mk(args...) = Expr(:call, args...)

@testset "XSPEC model string — Issue #187" begin

    # ── 1. XSPEC_MODEL_NAMES lookup table ────────────────────────────────────
    @testset "XSPEC_MODEL_NAMES table" begin
        # Every model that actually exists in this package must be present
        # Additive
        @test haskey(XSPEC_MODEL_NAMES, "powerlaw")
        @test haskey(XSPEC_MODEL_NAMES, "po")            # alias
        @test haskey(XSPEC_MODEL_NAMES, "cutoffpl")
        @test haskey(XSPEC_MODEL_NAMES, "bbody")
        @test haskey(XSPEC_MODEL_NAMES, "bb")            # alias
        @test haskey(XSPEC_MODEL_NAMES, "bremss")
        @test haskey(XSPEC_MODEL_NAMES, "kerrdisk")
        @test haskey(XSPEC_MODEL_NAMES, "kyrline")
        @test haskey(XSPEC_MODEL_NAMES, "laor")
        @test haskey(XSPEC_MODEL_NAMES, "diskline")
        @test haskey(XSPEC_MODEL_NAMES, "gaussian")
        @test haskey(XSPEC_MODEL_NAMES, "gaus")          # alias
        @test haskey(XSPEC_MODEL_NAMES, "jet")
        @test haskey(XSPEC_MODEL_NAMES, "optxagnf")
        @test haskey(XSPEC_MODEL_NAMES, "diskbb")
        # Multiplicative
        @test haskey(XSPEC_MODEL_NAMES, "phabs")
        @test haskey(XSPEC_MODEL_NAMES, "wndabs")
        @test haskey(XSPEC_MODEL_NAMES, "tbabs")
        # Convolutional
        @test haskey(XSPEC_MODEL_NAMES, "cflux")
        @test haskey(XSPEC_MODEL_NAMES, "kerrconv")

        # Aliases must resolve to the CORRECT Julia type
        @test XSPEC_MODEL_NAMES["po"]    === :XS_PowerLaw
        @test XSPEC_MODEL_NAMES["gaus"]  === :XS_Gaussian
        @test XSPEC_MODEL_NAMES["bb"]    === :XS_BlackBody

        # Correct Julia type names (the critical ones that differ from XSPEC names)
        @test XSPEC_MODEL_NAMES["tbabs"]   === :XS_NeutralHydrogenAbsorption
        @test XSPEC_MODEL_NAMES["phabs"]   === :XS_PhotoelectricAbsorption
        @test XSPEC_MODEL_NAMES["wndabs"]  === :XS_WarmAbsorption
        @test XSPEC_MODEL_NAMES["diskbb"]  === :XS_DiskBlackBody
        @test XSPEC_MODEL_NAMES["bremss"]  === :XS_BremsStrahlung
        @test XSPEC_MODEL_NAMES["cflux"]   === :XS_CalculateFlux

        # All values must be Symbols, all keys non-empty strings
        @test all(v isa Symbol  for v in values(XSPEC_MODEL_NAMES))
        @test all(!isempty(k)   for k in keys(XSPEC_MODEL_NAMES))

        # Models that do NOT exist in this package must NOT appear
        @test !haskey(XSPEC_MODEL_NAMES, "apec")    # not in XSPECModels.jl
        @test !haskey(XSPEC_MODEL_NAMES, "mekal")   # not in XSPECModels.jl
    end

    # ── 2. Single atomic model parsing ───────────────────────────────────────
    @testset "Single atomic models" begin
        @test parse_xspec_model_string("powerlaw") == _mk(:XS_PowerLaw)
        @test parse_xspec_model_string("bbody")    == _mk(:XS_BlackBody)
        @test parse_xspec_model_string("bremss")   == _mk(:XS_BremsStrahlung)
        @test parse_xspec_model_string("laor")     == _mk(:XS_Laor)
        @test parse_xspec_model_string("diskline") == _mk(:XS_DiskLine)
        @test parse_xspec_model_string("gaussian") == _mk(:XS_Gaussian)
        @test parse_xspec_model_string("diskbb")   == _mk(:XS_DiskBlackBody)
        @test parse_xspec_model_string("cutoffpl") == _mk(:XS_CutOffPowerLaw)
        @test parse_xspec_model_string("phabs")    == _mk(:XS_PhotoelectricAbsorption)
        @test parse_xspec_model_string("wndabs")   == _mk(:XS_WarmAbsorption)
        @test parse_xspec_model_string("tbabs")    == _mk(:XS_NeutralHydrogenAbsorption)
        @test parse_xspec_model_string("cflux")    == _mk(:XS_CalculateFlux)
        @test parse_xspec_model_string("kerrconv") == _mk(:XS_Kerrconv)
    end

    # ── 3. Aliases ────────────────────────────────────────────────────────────
    @testset "Aliases resolve to same type" begin
        @test parse_xspec_model_string("po")   == parse_xspec_model_string("powerlaw")
        @test parse_xspec_model_string("gaus") == parse_xspec_model_string("gaussian")
        @test parse_xspec_model_string("bb")   == parse_xspec_model_string("bbody")
    end

    # ── 4. Case insensitivity ─────────────────────────────────────────────────
    @testset "Case insensitivity" begin
        @test parse_xspec_model_string("POWERLAW")       == parse_xspec_model_string("powerlaw")
        @test parse_xspec_model_string("Phabs")          == parse_xspec_model_string("phabs")
        @test parse_xspec_model_string("PHABS*POWERLAW") == parse_xspec_model_string("phabs*powerlaw")
        @test parse_xspec_model_string("Tbabs*PowerLaw") == parse_xspec_model_string("tbabs*powerlaw")
    end

    # ── 5. Whitespace tolerance ───────────────────────────────────────────────
    @testset "Whitespace tolerance" begin
        @test parse_xspec_model_string("  phabs * powerlaw  ") ==
              parse_xspec_model_string("phabs*powerlaw")
        @test parse_xspec_model_string("phabs * ( powerlaw + bbody )") ==
              parse_xspec_model_string("phabs*(powerlaw+bbody)")
    end

    # ── 6. Binary multiplication ──────────────────────────────────────────────
    @testset "Multiplication M*A" begin
        result   = parse_xspec_model_string("phabs*powerlaw")
        expected = _mk(:*, _mk(:XS_PhotoelectricAbsorption), _mk(:XS_PowerLaw))
        @test result == expected

        result2   = parse_xspec_model_string("tbabs*powerlaw")
        expected2 = _mk(:*, _mk(:XS_NeutralHydrogenAbsorption), _mk(:XS_PowerLaw))
        @test result2 == expected2
    end

    # ── 7. Binary addition ────────────────────────────────────────────────────
    @testset "Addition A+A" begin
        result   = parse_xspec_model_string("powerlaw+bbody")
        expected = _mk(:+, _mk(:XS_PowerLaw), _mk(:XS_BlackBody))
        @test result == expected
    end

    # ── 8. CRITICAL: Operator precedence (* > +) ──────────────────────────────
    @testset "Operator precedence (* binds tighter than +)" begin
        # "phabs*powerlaw+bbody" must parse as (phabs*powerlaw)+bbody
        # NOT as phabs*(powerlaw+bbody)
        result   = parse_xspec_model_string("phabs*powerlaw+bbody")
        expected = _mk(:+,
            _mk(:*, _mk(:XS_PhotoelectricAbsorption), _mk(:XS_PowerLaw)),
            _mk(:XS_BlackBody))
        @test result == expected

        # Verify these are different
        without_parens = parse_xspec_model_string("phabs*powerlaw+bbody")
        with_parens    = parse_xspec_model_string("phabs*(powerlaw+bbody)")
        @test without_parens != with_parens   # MUST be different — astrophysically critical!
    end

    # ── 9. Parentheses override precedence ────────────────────────────────────
    @testset "Parentheses" begin
        result   = parse_xspec_model_string("phabs*(powerlaw+bbody)")
        expected = _mk(:*,
            _mk(:XS_PhotoelectricAbsorption),
            _mk(:+, _mk(:XS_PowerLaw), _mk(:XS_BlackBody)))
        @test result == expected
    end

    # ── 10. Left-associativity ────────────────────────────────────────────────
    @testset "Left-associativity" begin
        # powerlaw+bbody+gaussian → (powerlaw+bbody)+gaussian
        result   = parse_xspec_model_string("powerlaw+bbody+gaussian")
        expected = _mk(:+,
            _mk(:+, _mk(:XS_PowerLaw), _mk(:XS_BlackBody)),
            _mk(:XS_Gaussian))
        @test result == expected

        # phabs*tbabs*powerlaw → (phabs*tbabs)*powerlaw
        result2   = parse_xspec_model_string("phabs*tbabs*powerlaw")
        expected2 = _mk(:*,
            _mk(:*, _mk(:XS_PhotoelectricAbsorption), _mk(:XS_NeutralHydrogenAbsorption)),
            _mk(:XS_PowerLaw))
        @test result2 == expected2
    end

    # ── 11. Complex nesting ───────────────────────────────────────────────────
    @testset "Complex nested expressions" begin
        # phabs*(gaussian*(powerlaw+bbody))
        result = parse_xspec_model_string("phabs*(gaussian*(powerlaw+bbody))")
        @test result isa Expr
        @test result.head == :call && result.args[1] == :*
        @test result.args[2] == _mk(:XS_PhotoelectricAbsorption)
        # right subtree must be gaussian*(...)
        right = result.args[3]
        @test right.args[1] == :* && right.args[2] == _mk(:XS_Gaussian)

        # Deeply nested — must not error
        @test parse_xspec_model_string(
            "phabs*(wndabs*(powerlaw+bbody)+laor)"
        ) isa Expr
    end

    # ── 12. Real-world astrophysics model strings (all in this package) ───────
    @testset "Real-world model strings" begin
        # These are actual model strings used in X-ray astronomy
        real_models = [
            # classic absorbed power-law (AGN baseline)
            "phabs*powerlaw",
            # tbabs is NeutralHydrogenAbsorption in this package
            "tbabs*powerlaw",
            # double absorption column
            "phabs*tbabs*powerlaw",
            # soft excess
            "phabs*(powerlaw+bbody)",
            # disk + iron line
            "phabs*(diskbb+laor)",
            # multi-component
            "phabs*(powerlaw+bbody+gaussian)",
            # cutoff power-law
            "phabs*cutoffpl",
            # warm absorber
            "wndabs*powerlaw",
            # convolution flux calculation
            "cflux*powerlaw",
            # kerr convolution on disk
            "kerrconv*diskline",
            # bremsstrahlung with absorption
            "phabs*bremss",
            # multiple additive with different absorptions
            "phabs*(bbody+powerlaw+diskline)",
            # jet model
            "phabs*jet",
        ]
        for s in real_models
            result = parse_xspec_model_string(s)
            @test result isa Expr "Failed to parse: $s"
            @test result.head == :call "Wrong head for: $s"
        end
    end

    # ── 13. Convolutional model strings ───────────────────────────────────────
    @testset "Convolutional models" begin
        @test parse_xspec_model_string("cflux*powerlaw") ==
              _mk(:*, _mk(:XS_CalculateFlux), _mk(:XS_PowerLaw))
        @test parse_xspec_model_string("kerrconv*diskline") ==
              _mk(:*, _mk(:XS_Kerrconv), _mk(:XS_DiskLine))
    end

    # ── 14. Error handling ────────────────────────────────────────────────────
    @testset "Error handling" begin
        # Unknown model names
        @test_throws Exception parse_xspec_model_string("apec")      # not in package
        @test_throws Exception parse_xspec_model_string("mekal")     # not in package
        @test_throws Exception parse_xspec_model_string("notamodel")

        # Empty string
        @test_throws Exception parse_xspec_model_string("")
        @test_throws Exception parse_xspec_model_string("   ")

        # Structural parse errors
        @test_throws Exception parse_xspec_model_string("*powerlaw")        # leading *
        @test_throws Exception parse_xspec_model_string("powerlaw*")        # trailing *
        @test_throws Exception parse_xspec_model_string("+powerlaw")        # leading +
        @test_throws Exception parse_xspec_model_string("phabs*(powerlaw")  # missing )
        @test_throws Exception parse_xspec_model_string("phabs**powerlaw")  # double *

        # Error message must mention the offending name
        err_msg = try
            parse_xspec_model_string("xyz_not_a_model")
            ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("xyz_not_a_model", err_msg)
    end

    # ── 15. Reverse table _XSPEC_REVERSE consistency ──────────────────────────
    @testset "Reverse lookup table consistency" begin
        # Every entry in _XSPEC_REVERSE must point back to a valid forward entry
        for (jsym, xname) in XSPECModels._XSPEC_REVERSE
            @test haskey(XSPEC_MODEL_NAMES, xname)
            @test XSPEC_MODEL_NAMES[xname] === jsym
        end
        # Key types are present
        @test haskey(XSPECModels._XSPEC_REVERSE, :XS_PowerLaw)
        @test haskey(XSPECModels._XSPEC_REVERSE, :XS_PhotoelectricAbsorption)
        @test haskey(XSPECModels._XSPEC_REVERSE, :XS_NeutralHydrogenAbsorption)
        @test haskey(XSPECModels._XSPEC_REVERSE, :XS_CalculateFlux)

        # Canonical form must be the longer one (not alias)
        @test XSPECModels._XSPEC_REVERSE[:XS_PowerLaw] == "powerlaw"    # not "po"
        @test XSPECModels._XSPEC_REVERSE[:XS_Gaussian] == "gaussian"    # not "gaus"
        @test XSPECModels._XSPEC_REVERSE[:XS_BlackBody] == "bbody"      # not "bb"
    end

    # ── 16. xspec_model_string — reverse direction ────────────────────────────
    @testset "xspec_model_string (reverse)" begin
        # Atomic models
        @test xspec_model_string(XS_PowerLaw())                 == "powerlaw"
        @test xspec_model_string(XS_BlackBody())                == "bbody"
        @test xspec_model_string(XS_BremsStrahlung())           == "bremss"
        @test xspec_model_string(XS_Laor())                     == "laor"
        @test xspec_model_string(XS_DiskLine())                 == "diskline"
        @test xspec_model_string(XS_Gaussian())                 == "gaussian"
        @test xspec_model_string(XS_DiskBlackBody())            == "diskbb"
        @test xspec_model_string(XS_CutOffPowerLaw())           == "cutoffpl"  # canonical, not "zcutoffpl"
        @test xspec_model_string(XS_PhotoelectricAbsorption())  == "phabs"
        @test xspec_model_string(XS_WarmAbsorption())           == "wndabs"
        @test xspec_model_string(XS_NeutralHydrogenAbsorption()) == "tbabs"
        @test xspec_model_string(XS_CalculateFlux())            == "cflux"
        @test xspec_model_string(XS_Kerrconv())                 == "kerrconv"

        # Composite: multiplication
        m_mul = XS_PhotoelectricAbsorption() * XS_PowerLaw()
        @test xspec_model_string(m_mul) == "phabs*powerlaw"

        # Composite: addition
        m_add = XS_PowerLaw() + XS_BlackBody()
        @test xspec_model_string(m_add) == "powerlaw+bbody"

        # Composite: multiplication with addition (parens needed!)
        m_complex = XS_PhotoelectricAbsorption() * (XS_PowerLaw() + XS_BlackBody())
        @test xspec_model_string(m_complex) == "phabs*(powerlaw+bbody)"

        # Chained multiplication (no parens needed)
        m_chain = XS_PhotoelectricAbsorption() * XS_NeutralHydrogenAbsorption() * XS_PowerLaw()
        @test xspec_model_string(m_chain) == "phabs*tbabs*powerlaw"

        # Three-way addition
        m_three = XS_PowerLaw() + XS_BlackBody() + XS_Gaussian()
        @test xspec_model_string(m_three) == "powerlaw+bbody+gaussian"
    end

    # ── 17. Round-trip: parse → eval → xspec_model_string → parse ─────────────
    @testset "Round-trip: string → Expr → model → string → Expr" begin
        roundtrip_cases = [
            "powerlaw",
            "phabs*powerlaw",
            "tbabs*powerlaw",
            "powerlaw+bbody",
            "phabs*(powerlaw+bbody)",
            "phabs*tbabs*powerlaw",
            "phabs*(powerlaw+bbody+gaussian)",
            "cflux*powerlaw",
            "kerrconv*diskline",
        ]
        for s in roundtrip_cases
            # Parse to model
            model = eval(parse_xspec_model_string(s))
            # Convert back to string
            s2    = xspec_model_string(model)
            # Re-parse and check AST equality
            expr1 = parse_xspec_model_string(s)
            expr2 = parse_xspec_model_string(s2)
            @test expr1 == expr2   "Round-trip failed for: $s → $s2"
        end
    end

    # ── 18. @xspec_str macro compiles to the right types ─────────────────────
    @testset "@xspec_str macro" begin
        # Single model
        m1 = xspec"powerlaw"
        @test m1 isa XS_PowerLaw

        # Multiplication
        m2 = xspec"phabs*powerlaw"
        @test m2 isa CompositeModel
        @test m2.left  isa XS_PhotoelectricAbsorption
        @test m2.right isa XS_PowerLaw

        # Absorption on whole sum
        m3 = xspec"tbabs*(powerlaw+bbody)"
        @test m3 isa CompositeModel
        @test m3.left  isa XS_NeutralHydrogenAbsorption
        @test m3.right isa CompositeModel
        @test m3.right.left  isa XS_PowerLaw
        @test m3.right.right isa XS_BlackBody

        # Precedence: phabs*powerlaw+bbody → (phabs*powerlaw)+bbody
        m4 = xspec"phabs*powerlaw+bbody"
        @test m4 isa CompositeModel
        @test m4.left  isa CompositeModel   # phabs*powerlaw is the LEFT
        @test m4.right isa XS_BlackBody     # bbody is the RIGHT

        # cflux convolution
        m5 = xspec"cflux*powerlaw"
        @test m5 isa CompositeModel
        @test m5.left isa XS_CalculateFlux
        @test m5.right isa XS_PowerLaw

        # Case insensitive
        m6 = xspec"PHABS*PowerLaw"
        @test m6 isa CompositeModel
        @test m6.left  isa XS_PhotoelectricAbsorption
        @test m6.right isa XS_PowerLaw
    end

end   # @testset "XSPEC model string — Issue #187"

println("\n✅  All XSPEC model string tests passed!")
