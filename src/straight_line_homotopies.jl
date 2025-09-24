export straight_line_homotopy,
    linear_path,
    specified_system


# constructing a straight line homotopy from G (when t=0) to F (when t=1).
function straight_line_homotopy(F, G, t)
    n = length(F)
    HR = parent(t)
    gamma = CCi(rand(Float64), rand(Float64))
    H = zeros(HR, 0)
    for i in 1:n
        H = push!(H, (1 - t) * gamma * G[i] + t * F[i])
    end
    H
end

function max_degree(H::Matrix)
    generic_F = evaluate_matrix(H, CCi(rand(ComplexF64)))
#    return map(i -> degrees(generic_F[i]), 1:length(generic_F))
    return map(i -> maximum(map(j -> sum(degrees(j)), AbstractAlgebra.monomials(generic_F[i]))), 1:length(generic_F))
end

# constructing a linear path from p0 to p1
function linear_path(p0, p1, t)
    n = length(p0)
    HR = parent(t)
    p = zeros(HR, n)
    for i = 1:n
        p[i] = p0[i] * (1 - t) + p1[i] * t
    end
    p
end

# constructing a specified_system from p0 to p1
function specified_system(p0, p1, F)
    HR = base_ring(F[1])
    n = length(F)
    t = gens(HR)[1]
    p = linear_path(p0, p1, t)
#    if length(p) == 1
#        p = p[1]
#    end
    Fp = zeros(HR, n)
    for i = 1:n
        Fp[i] = AbstractAlgebra.evaluate(F[i], p)
    end
    hcat(Fp)
end
