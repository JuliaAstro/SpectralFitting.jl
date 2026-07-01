using SpectralFitting, Test


dummy_data1 = make_dummy_dataset((E) -> E^(-3.0) + E^(-2.3); units = u"counts / (s * keV)")

model1 = PowerLaw() + PowerLaw()
model2 = sum(PowerLaw() for i in 1:3)

prob = FittingProblem(model1 => dummy_data1, model2 => dummy_data1)

mapping = SpectralFitting._build_parameter_mapping(prob)
@test mapping == ([1, 2, 3, 4], [5, 6, 7, 8, 9, 10])

# now we bind a parameter
bind!(prob, (1, :a1, :K) => (2, :a1, :K))

mapping = SpectralFitting._build_parameter_mapping(prob)
@test mapping == ([1, 2, 3, 4], [1, 6, 7, 8, 9, 10])

model1 = PowerLaw()
model2 = PowerLaw() + PowerLaw()

# binding weird model combinations
prob = FittingProblem(model1 => dummy_data1, model2 => dummy_data1)
bind!(prob, (1, :K) => (2, :a2, :K))

mapping = SpectralFitting._build_parameter_mapping(prob)
@test mapping == ([1, 2], [3, 4, 1, 6])

model1 =
    PowerLaw(K = FitParam(3.0), a = FitParam(4.0)) +
    PowerLaw(K = FitParam(1.0), a = FitParam(2.0))
model2 =
    PowerLaw(K = FitParam(9.0), a = FitParam(10.0)) +
    PowerLaw(K = FitParam(7.0), a = FitParam(8.0)) +
    PowerLaw(K = FitParam(5.0), a = FitParam(6.0))

prob = FittingProblem(model1 => dummy_data1, model2 => dummy_data1)
bind!(prob, (1, :a1, :K) => (2, :a1, :K))
bind!(prob, (1, :a2, :K) => (2, :a2, :K))

mapping = SpectralFitting._build_parameter_mapping(prob)
@test mapping == ([1, 2, 3, 4], [1, 6, 3, 8, 9, 10])

model1 = PowerLaw()
model2 = PowerLaw() + PowerLaw()

# test multiple bindings with multiple models and multiple datasets
prob = FittingProblem(model1 => dummy_data1, model1 => dummy_data1, model1 => dummy_data1)
bind!(prob, (1, :K) => (2, :K) => (3, :K))
bind!(prob, (1, :a) => (2, :a) => (3, :a))
mapping = SpectralFitting._build_parameter_mapping(prob)
@test mapping == ([1, 2], [1, 2], [1, 2])

prob = FittingProblem(model1 => dummy_data1, model1 => dummy_data1, model1 => dummy_data1)
bindall!(prob, :K)
bindall!(prob, :a)
mapping = SpectralFitting._build_parameter_mapping(prob)
@test mapping == ([1, 2], [1, 2], [1, 2])

# 3 models, 1 bound parameter
prob = FittingProblem(model1 => dummy_data1, model1 => dummy_data1, model1 => dummy_data1)
bindall!(prob, :a)
mapping = SpectralFitting._build_parameter_mapping(prob)
@test mapping == ([1, 2], [3, 2], [5, 2])

# 3 models (2 the same, 1 different), bind one parameter between the two
prob = FittingProblem(model1 => dummy_data1, model2 => dummy_data1, model1 => dummy_data1)
bind!(prob, (1, :a) => (2, :a2, :a) => (3, :a))
mapping = SpectralFitting._build_parameter_mapping(prob)
@test mapping == ([1, 2], [3, 4, 5, 2], [7, 2])

# 3 models, 1 frozen parameter (that needs to be skipped), 1 bound parameter
model1 = PowerLaw() + PowerLaw()
model1.a1.K.frozen = true
prob = FittingProblem(model1 => dummy_data1, model1 => dummy_data1, model1 => dummy_data1)
bindall!(prob, (:a1, :a), (:a2, :a))

# check that the free parameter binding adjustment works okay
mapping = SpectralFitting._build_parameter_mapping(prob)
@test mapping == ([1, 2, 3, 4], [5, 2, 7, 4], [9, 2, 11, 4])

# does this construct the correct parameter cache

# 3 models, 1 frozen parameter (that doesn't need to be skipped), 1 bound parameter
model1.a1.K.frozen = false
model1.a2.a.frozen = true
prob = FittingProblem(model1 => dummy_data1, model1 => dummy_data1, model1 => dummy_data1)
bindall!(prob, (:a1, :K), (:a1, :a))


mapping = SpectralFitting._build_parameter_mapping(prob)
@test mapping == ([1, 2, 3, 4], [1, 2, 7, 8], [1, 2, 11, 12])

# can we bind parameters within the same model
prob = FittingProblem(model1 => dummy_data1)
bind!(prob, (1, :a1, :K) => (1, :a1, :a))
mapping = SpectralFitting._build_parameter_mapping(prob)
@test mapping == ([1, 1, 3, 4],)

# can we bind parameters within the same model
model1.a2.a.frozen = false
prob = FittingProblem(model1 => dummy_data1, model1 => dummy_data1)
bind!(prob, (1, :a1, :K) => (1, :a1, :a))
bind!(prob, (1, :a1, :K) => (2, :a1, :K))

bind!(prob, (1, :a1, :a) => (2, :a1, :a))

bind!(prob, (1, :a2, :a) => (2, :a2, :K))
SpectralFitting.simplify!(prob)
mapping = SpectralFitting._build_parameter_mapping(prob)
@test mapping == ([1, 1, 3, 4], [1, 1, 4, 8])
