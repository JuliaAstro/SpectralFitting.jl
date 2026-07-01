using Test
using SpectralFitting


struct EvalCountingModel{D, T} <: AbstractTableModel{T, Additive}
    table::D
    K::T
    a::T
end

function EvalCountingModel(; K = FitParam(1.0), a = FitParam(3.0))
    return EvalCountingModel(Int[0], K, a)
end

function SpectralFitting.invoke!(output, domain, model::EvalCountingModel)
    model.table[1] += 1
    return @. output = domain[1:(end - 1)] .+ model.a
end

domain = collect(range(0.0, 10.0, 100))

model = AutoCache(EvalCountingModel())

# running the model several times should only hit the counter once
for i in 1:100
    invokemodel(domain, model)
end
@test model.model.table[1] == 1

# changing the parameter should hit the counter again
set_value!(model.model.a, 5.0)
for i in 1:100
    invokemodel(domain, model)
end
@test model.model.table[1] == 2

# modifying the domain should change the cache as well
domain = collect(range(0.1, 5.0, 10))
for i in 1:100
    invokemodel(domain, model)
end
@test model.model.table[1] == 3

# now as a composite
domain = collect(range(0.0, 10.0, 100))
model = DummyMultiplicative() * AutoCache(EvalCountingModel())

# running the model several times should only hit the counter once
for i in 1:10
    invokemodel(domain, model)
end
@test model.a1.model.table[1] == 1

# changing the parameter should hit the counter again
set_value!(model.a1.a, 5.0)
for i in 1:100
    invokemodel(domain, model)
end
@test model.a1.model.table[1] == 2

# but not if it's normalisation
set_value!(model.a1.K, 10.0)
for i in 1:100
    invokemodel(domain, model)
end
@test model.a1.model.table[1] == 2

# modifying the domain should change the cache as well
domain = collect(range(0.1, 5.0, 10))
for i in 1:100
    invokemodel(domain, model)
end
@test model.a1.model.table[1] == 3
