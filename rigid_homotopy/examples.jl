include("rigid_homotopy.jl")

# a two-variable example
CCi = AcbField()
eR, (x,y,η) = CCi["x","y","η"]
HR, (t) = eR["t"]

F = [x^2 + x*y+ 2*y^2 - 4  x^2 + y^2 - 1]
point = [CCi(-.576287,-1.04024), CCi(-1.39128,.430883)]

M = rand(ComplexF64, 2, 2)
A1 = matrix(convert_to_box_matrix((M- M')/2,CCi))
M = rand(ComplexF64, 2, 2)
A2 = matrix(convert_to_box_matrix((M- M')/2,CCi))

w1 = matrix_exponential(A1, t; order = 5)
w2 = matrix_exponential(A2, t; order = 5)
wt = [A1, A2]
p, iters = rigid_track(F, point, .1, wt)


# comparison with the certified linear homotopy
# These are needed to construct the linear homotopy
R = parent(F[1])
n = length(F)
vars_F = gens(R)[1:n]

w_t = map(i -> evaluate_matrix(Matrix(i), 1), wt)
transformations = map(i -> i*vars_F, w_t)
target_system = map(i -> AbstractAlgebra.evaluate(F[i], [transformations[i]; 0]), 1:n)

H = Matrix((1-t)*matrix(transpose(F))+ t*matrix(target_system))
HH = Matrix((1-t)*matrix(F)+ t*matrix(transpose(target_system)))
tracking_without_predictor(HH, point, .1; iterations_count = true) # 55 iterations
track(H, point, .1; iterations_count = true) # 32 iterations




# a two-variable example
CCi = AcbField()
eR, (x,y,η) = CCi["x","y","η"]
HR, (t) = eR["t"]

A_0 = matrix(convert_to_box_matrix(rand(ComplexF64, 2, 2),CCi))
A_1 = matrix(convert_to_box_matrix(rand(ComplexF64, 2, 2),CCi))
A_2 = matrix(convert_to_box_matrix(rand(ComplexF64, 2, 2),CCi))
B_0 = matrix(convert_to_box_matrix(rand(ComplexF64, 2, 2),CCi))
B_1 = matrix(convert_to_box_matrix(rand(ComplexF64, 2, 2),CCi))
B_2 = matrix(convert_to_box_matrix(rand(ComplexF64, 2, 2),CCi))


F = [CCi(.0665973,+.715566)*x^2+CCi(-.443541,+.0448575)*x*y+CCi(.320556,-.232718)*y^2+CCi(.854157,-1.0638)*x+CCi(.827374,+.264383)*y-CCi(.494373,.295369)
    CCi(-1.10874,-.0359951)*x^2+CCi(-1.12079,-.558265)*x*y+CCi(-.255562,-1.11227)*y^2+CCi(1.59943,-.472601)*x+CCi(.327772,+.644435)*y-CCi(.330381,+.00997843)]
F = Matrix(transpose(F))
points = [[CCi(1.07906,-.726546), CCi(.306091,+1.1297)], [CCi(.254121,+1.04198), CCi(-1.39188,+.148427)],
      [CCi(1.00155,.646277), CCi(.278009,-.790194)], [CCi(.265832,+.0274612), CCi(.288628,-.168463)]]

M = rand(ComplexF64, 2, 2)
A1 = matrix(convert_to_box_matrix((M- M')/2,CCi))
M = rand(ComplexF64, 2, 2)
A2 = matrix(convert_to_box_matrix((M- M')/2,CCi))

w1 = matrix_exponential(A1, t; order = 1)
w2 = matrix_exponential(A2, t; order = 1)
At = [A1, A2]
wt = [w1, w2]
p1, iters1 = rigid_track2(F, points[1], .1, wt, At)
p2, iters2 = rigid_track(F, points[2], .1, wt)
p3, iters3 = rigid_track(F, points[3], .1, wt)
p4, iters4 = rigid_track(F, points[4], .1, wt)
iters1 + iters2 + iters3 + iters4

# comparison with the certified linear homotopy
# These are needed to construct the linear homotopy
    R = parent(F[1])
n = length(F)
vars_F = gens(R)[1:n]

w_t = map(i -> evaluate_matrix(Matrix(i), 1), wt)
transformations = map(i -> i*vars_F, w_t)
target_system = map(i -> AbstractAlgebra.evaluate(F[i], [transformations[i]; 0]), 1:n)

H = Matrix((1-t)*matrix(transpose(F))+ t*matrix(target_system))
HH = Matrix((1-t)*CCi(exp(im*rand(Int)))*matrix(F)+ t*matrix(transpose(target_system)))
lp1, liter1 = tracking_without_predictor(HH, points[1], .1; iterations_count = true)
lp2, liter2 = tracking_without_predictor(HH, points[2], .1; iterations_count = true)
lp3, liter3 = tracking_without_predictor(HH, points[3], .1; iterations_count = true)
lp4, liter4 = tracking_without_predictor(HH, points[4], .1; iterations_count = true)
liter1 + liter2 + liter3 + liter4
track(H, point, .1; iterations_count = true) # 32 iterations






# an example with 4 variables
CCi = AcbField()
eR, (x,y,z,w,η) = CCi["x","y","z","w","η"]
HR, (t) = eR["t"]

# Define system F with 4 equations
F = [x^2+x*y+y*z+2*w^2-4 x^2+y^2+z*w-1 z^2+w*x+y*z-2 w^2+x*z+y*w+3]
point = [CCi(1.73292,+.199961), CCi(.466897,-1.46505), CCi(-1.27595,+.55608), CCi(.216577,-.43465)]

# Create random skew-Hermitian matrices for unitary paths
M1 = rand(ComplexF64, 4, 4)
M2 = rand(ComplexF64, 4, 4)
M3 = rand(ComplexF64, 4, 4)
M4 = rand(ComplexF64, 4, 4)

# Convert to skew-Hermitian matrices (A = -A*)
A1 = matrix(convert_to_box_matrix((M1 - M1')/2, CCi))
A2 = matrix(convert_to_box_matrix((M2 - M2')/2, CCi))
A3 = matrix(convert_to_box_matrix((M3 - M3')/2, CCi))
A4 = matrix(convert_to_box_matrix((M4 - M4')/2, CCi))

# Create unitary matrices using matrix exponential
w1 = matrix_exponential(A1, t; order = 10)
w2 = matrix_exponential(A2, t; order = 10)
w3 = matrix_exponential(A3, t; order = 10)
w4 = matrix_exponential(A4, t; order = 10)

# Collect all unitary matrices
wt = [w1, w2, w3, w4]


p, iters = rigid_track(F, point, .1, wt)
iters # 11 iterations


# comparison with the certified linear homotopy
# These are needed to construct the linear homotopy
R = parent(F[1])
n = length(F)
vars_F = gens(R)[1:n]

w_t = map(i -> evaluate_matrix(Matrix(i), 1), wt)
transformations = map(i -> i*vars_F, w_t)
target_system = map(i -> AbstractAlgebra.evaluate(F[i], [transformations[i]; 0]), 1:n)

H = Matrix((1-t)*matrix(transpose(F))+ t*matrix(target_system))
HH = Matrix((1-t)*matrix(F)+ t*matrix(transpose(target_system)))
tracking_without_predictor(HH, point, .1; iterations_count = true) # 55 iterations
track(H, point, .1; iterations_count = true) # 32 iterations


# an example with 4 variables of higher degree
CCi = AcbField()
eR, (x,y,z,w,η) = CCi["x","y","z","w","η"]
HR, (t) = eR["t"]

# Define system F with 4 equations
F = [ x^4 + 3*x^2*y^2 - 2*y^3*z + z^3*w + x*y*z*w - y^2*w^2 + x^2*z^2 - 5,
    y^4 - x^2*z^3 + 4*w^3 - z^4 + y*z^2*w + x*y*w^2 - x^3*w + 3,
    z^4*w + y^3*z^2 + x^2*y*z + w^4 - x*z^2*w + x*y*z^2 - z*w^2 + 1,
    x^3*y^2 - z^3*w + y^2*z^3 + x*w^2 + x^2*y^2 - x*y*z*w + y*w^3 - 6]
F = Matrix(transpose(F))
point = [CCi(-1.69158,+.583695), CCi(.931417,+.600535), CCi(.38314,+1.48638), CCi(.894437,-.89492)]

# Create random skew-Hermitian matrices for unitary paths
M1 = rand(ComplexF64, 4, 4)
M2 = rand(ComplexF64, 4, 4)
M3 = rand(ComplexF64, 4, 4)
M4 = rand(ComplexF64, 4, 4)

# Convert to skew-Hermitian matrices (A = -A*)
A1 = matrix(convert_to_box_matrix((M1 - M1')/2, CCi))
A2 = matrix(convert_to_box_matrix((M2 - M2')/2, CCi))
A3 = matrix(convert_to_box_matrix((M3 - M3')/2, CCi))
A4 = matrix(convert_to_box_matrix((M4 - M4')/2, CCi))

# Create unitary matrices using matrix exponential
w1 = matrix_exponential(A1, t; order = 10)
w2 = matrix_exponential(A2, t; order = 10)
w3 = matrix_exponential(A3, t; order = 10)
w4 = matrix_exponential(A4, t; order = 10)

# Collect all unitary matrices
wt = [w1, w2, w3, w4]


p, iters= rigid_track(F, point, .1, wt)
iters # 11 iterations


# comparison with the certified linear homotopy
# These are needed to construct the linear homotopy
R = parent(F[1])
n = length(F)
vars_F = gens(R)[1:n]

w_t = map(i -> evaluate_matrix(Matrix(i), 1), wt)
transformations = map(i -> i*vars_F, w_t)
target_system = map(i -> AbstractAlgebra.evaluate(F[i], [transformations[i]; 0]), 1:n)

H = Matrix((1-t)*matrix(transpose(F))+ t*matrix(target_system))
HH = Matrix((1-t)*matrix(F)+ t*matrix(transpose(target_system)))
tracking_without_predictor(HH, point, .1; iterations_count = true) # 55 iterations