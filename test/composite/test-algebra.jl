using Test
using SpectralFitting


add_model = DummyAdditive()
mul_model = DummyMultiplicative()

# allowed operations
cm = add_model + add_model
@test modelkind(cm) == SpectralFitting.Additive()
cm = mul_model * add_model
@test modelkind(cm) == SpectralFitting.Additive()
cm = mul_model * mul_model
@test modelkind(cm) == SpectralFitting.Multiplicative()

# don't care so much about the error message just the error
@test_throws "" add_model * add_model
@test_throws "" add_model * mul_model
@test_throws "" mul_model + mul_model
