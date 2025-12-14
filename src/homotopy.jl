export straight_line_homotopy, linear_path, specified_system, max_degree

# ------------------------------------------------------------------------------
# Straight Line Homotopy
# ------------------------------------------------------------------------------
function straight_line_homotopy(F, G, t)
    n = length(F)
    
    R = parent(F[1])
    CC = base_ring(R) 
    
    gamma = CC(rand(Float64), rand(Float64))
    
    HR = parent(t)
    H = zeros(HR, 0)
    for i in 1:n
        H = push!(H, (1 - t) * gamma * G[i] + t * F[i])
    end
    H
end

# ------------------------------------------------------------------------------
# Max Degree Calculation
# ------------------------------------------------------------------------------
function max_degree(H::Matrix)
    HR = parent(H[1])       # R[t]
    R  = base_ring(HR)      # R
    CC = base_ring(R)       # CC (AcbField)
    
    rand_t = CC(rand(ComplexF64))
    generic_F = evaluate_matrix(H, rand_t)

    return map(i -> maximum(map(j -> sum(degrees(j)), AbstractAlgebra.monomials(generic_F[i]))), 1:length(generic_F))
end

# ------------------------------------------------------------------------------
# Linear Path Construction
# ------------------------------------------------------------------------------
function linear_path(p0, p1, t)
    n = length(p0)
    HR = parent(t)
    p = zeros(HR, n)
    for i = 1:n
        p[i] = p0[i] * (1 - t) + p1[i] * t
    end
    p
end

# ------------------------------------------------------------------------------
# Specified System Construction
# ------------------------------------------------------------------------------
function specified_system(p0, p1, F)
    
    HR = base_ring(F[1]) 
    
    n = length(F)
    t = gens(HR)[1] 
    
    p = linear_path(p0, p1, t)

    Fp = zeros(HR, n)
    for i in 1:n
        Fp[i] = AbstractAlgebra.evaluate(F[i], p)
    end
    hcat(Fp)
end