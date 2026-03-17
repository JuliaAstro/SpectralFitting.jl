using Test
using SpectralFitting




function showstring(item)
    buffer = IOBuffer()
    Base.show(IOContext(buffer, :color => false), MIME"text/plain"(), item)
    strip(String(take!(buffer)))
end

# printing single
model = DummyAdditive()
string = showstring(model)
expected = """
┌ DummyAdditive
│      K ->  1 ± 0.1  ∈ [ 0, Inf ]   FREE
│      a ->  1 ± 0.1  ∈ [ 0, Inf ]   FREE
│      b ->  5                       FROZEN
└"""
@test string == expected

# printing composite
model = DummyMultiplicative() * (DummyAdditive() + DummyAdditive())
# test destructuring works
# test the output string is correct
string = showstring(model)
expected = """┌ CompositeModel with 3 model components:
│      m1 * (a2 + a1)
│ Model key and parameters:
│    a1 => DummyAdditive
│      K_1 ->  1 ± 0.1  ∈ [ 0, Inf ]   FREE
│      a_1 ->  1 ± 0.1  ∈ [ 0, Inf ]   FREE
│      b_1 ->  5                       FROZEN
│    a2 => DummyAdditive
│      K_2 ->  1 ± 0.1  ∈ [ 0, Inf ]   FREE
│      a_2 ->  1 ± 0.1  ∈ [ 0, Inf ]   FREE
│      b_2 ->  5                       FROZEN
│    m1 => DummyMultiplicative
│      a_3 ->  1 ± 0.1  ∈ [ 0, Inf ]   FREE
│      b_3 ->  5                       FROZEN
└"""
@test string == expected


# check string for more complex model
model =
    (DummyMultiplicative() * DummyMultiplicative() * DummyAdditive()) +
    DummyMultiplicative() * DummyAdditive()
string = showstring(model)
expected = """┌ CompositeModel with 5 model components:
│      (m3 * m2) * a2 + m1 * a1
│ Model key and parameters:
│    a1 => DummyAdditive
│      K_1 ->  1 ± 0.1  ∈ [ 0, Inf ]   FREE
│      a_1 ->  1 ± 0.1  ∈ [ 0, Inf ]   FREE
│      b_1 ->  5                       FROZEN
│    m1 => DummyMultiplicative
│      a_2 ->  1 ± 0.1  ∈ [ 0, Inf ]   FREE
│      b_2 ->  5                       FROZEN
│    a2 => DummyAdditive
│      K_2 ->  1 ± 0.1  ∈ [ 0, Inf ]   FREE
│      a_3 ->  1 ± 0.1  ∈ [ 0, Inf ]   FREE
│      b_3 ->  5                       FROZEN
│    m2 => DummyMultiplicative
│      a_4 ->  1 ± 0.1  ∈ [ 0, Inf ]   FREE
│      b_4 ->  5                       FROZEN
│    m3 => DummyMultiplicative
│      a_5 ->  1 ± 0.1  ∈ [ 0, Inf ]   FREE
│      b_5 ->  5                       FROZEN
└"""
@test string == expected

# test problem printing
model = DummyMultiplicative() * DummyAdditive()
data = make_dummy_dataset(x -> x)
prob = FittingProblem(model, data)
string = showstring(prob)
expected = """┌ FittingProblem:
│   . Models     : 1
│   . Datasets   : 1
│   Parameter Summary:
│   . Total      : 5
│   . Frozen     : 2
│   . Bound      : 0
│   . Free       : 3
└"""
@test string == expected
