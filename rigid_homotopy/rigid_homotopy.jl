include("../src/certified_monodromy_computation.jl")

export matrix_exponential, hermitian_conjugate, rigid_track

function matrix_exponential(A, t; order::Int = 10)
    n = size(A)[1]
    id_matrix = matrix(Matrix(1*I, n,n))

    res = id_matrix
    for i in 1:order
        res += (t^i / factorial(i)) * A^i
    end
    return res
end

function hermitian_conjugate(A, CCi)
    n = size(A)[1]
    res = convert_to_box_matrix(zeros(ComplexF64, n, n),CCi)
    for i in 1:n
        for j in 1:n
            ent = A[i,j]
            res[j,i] = CCi(real(ent), -imag(ent))
        end
    end
    return res
end


function rigid_track(F, point, r, w)

    R = parent(F[1])
    n = length(F)
    vars_F = gens(R)[1:n]
    h = 1/10
    tval = 0
    iter = 0

    while tval < 1
        iter += 1
        print("\rprogress t: $(round(tval, digits = 10))")
        print("\riter: $iter")
        w_t = map(i -> evaluate_matrix(Matrix(i), tval), w)
        transformations = map(i -> i*vars_F, w_t)
        transformed_F = map(i -> AbstractAlgebra.evaluate(F[i], [transformations[i]; 0]), 1:n)
        A = jacobian_inverse(matrix(transformed_F), point)

        point,r,A = refine_moore_box(Matrix(transpose(matrix(transformed_F))), point, r, A, 1/8)
        while krawczyk_test(Matrix(transpose(matrix(transformed_F))), point, r, A, 7/8) == false
            h = 1/2 * h
            midt = tval + h/2
            radii = h/2
            T = CCi("$midt +/- $radii")

            w_t = map(i -> evaluate_matrix(Matrix(i), T), w)
            transformations = map(i -> i*vars_F, w_t)
            transformed_F = map(i -> AbstractAlgebra.evaluate(F[i], [transformations[i]; 0]), 1:n)
        end

        tval = tval + h
    end

    w_t = map(i -> evaluate_matrix(Matrix(i), 1), w)
    transformations = map(i -> i*vars_F, w_t)
    transformed_F = map(i -> AbstractAlgebra.evaluate(F[i], [transformations[i]; 0]), 1:n)

    point,r,A = refine_moore_box(Matrix(transpose(matrix(transformed_F))), point, r, A, 1/100)
    return point,iter
end
