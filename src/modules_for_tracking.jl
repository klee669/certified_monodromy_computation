export krawczyk_operator_taylor_model,
       proceeding_step,
       refine_step,
       linear_tracking,
       hermite_tracking


# --------------------------------------------------------------------------
# Krawczyk operator in the Taylor model
function krawczyk_operator_taylor_model(
    H::Union{Matrix, Vector}, 
    lp::Vector, 
    tval::Number,
    A::AcbMatrix,
    r::Number,
)
    n   = length(H)
    eR  = base_ring(H[1])
    CC  = coefficient_ring(eR)
    η   = gens(eR)[end]
    lp  = push!(lp, η)

    eH    = evaluate_matrix(hcat(H), tval + η)
    ejac  = jac(eH)
    eHjac = zeros(eR, n, n)

    ball  = CC("+/- 1", "+/-1")
    mat   = [ball for _ in 1:n] |> x -> vcat(x, [CC(0)])

    for i in 1:n
        eH[i] = AbstractAlgebra.evaluate(eH[i], lp)
        for j in 1:n
            eHjac[i, j] = AbstractAlgebra.evaluate(ejac[i, j], lp + r * mat)
        end
    end

    I = identity_matrix(CC, n)
    (-1 / r) * A * eH + (I - A * matrix(eHjac)) * matrix(mat[1:n])
end


# --------------------------------------------------------------------------
# Shrinks step size `h` until the Krawczyk test passes
function proceeding_step(
    h::Number, 
    CCi::AcbField, 
    n::Int, 
    tm::AbstractAlgebra.Generic.MatSpaceElem, 
    K::AcbMatrix,
)
    while max_norm(K) > 7/8
        h /= 2
        radii = h / 2

        if h < 1e-100
            error("Step size h too small (underflow)")
        end

        T = CCi("$radii +/- $radii")
        input = zeros(CCi, n + 1)
        input[n + 1] = T

        K = evaluate_matrix(tm, input)
    end

    return h
end


# --------------------------------------------------------------------------
# Refine a point `x` at homotopy time `t`
function refine_step(
    H::Union{Matrix, Vector}, 
    t::Number, 
    x::Vector, 
    r::Number, 
    A::AcbMatrix, 
    h::Number; 
    threshold = 1 / 8,
)
    Ft = evaluate_matrix(Matrix(transpose(hcat(H))), t)
    x, r, A = refine_moore_box(Ft, x, r, A, threshold)

    v      = speed_vector(H, x, t, A)
    h     *= 5 / 4
    radii  = h / 2

    return Ft, x, r, A, v, h, radii
end


# --------------------------------------------------------------------------
# First step using linear predictor
function linear_tracking(
    H::Union{Matrix, Vector}, 
    t::Number, 
    x::Vector, 
    r::Number, 
    A::AcbMatrix, 
    h::Number, 
    n::Int, 
    refinement_threshold::Number,
)
    Ft, x, r, A, v, h, radii = refine_step(H, t, x, r, A, h; threshold = refinement_threshold)
    X = linear_predictor(H, v, x)

    tm = krawczyk_operator_taylor_model(H, X, t, A, r)

    T = CCi("$radii +/- $radii")
    input = zeros(CCi, n + 1)
    input[n + 1] = T

    K = evaluate_matrix(tm, input)
    h = proceeding_step(h, CCi, n, tm, K)

    return x, v, h, X, r, A
end


# --------------------------------------------------------------------------
# Step using Hermite predictor (with history)
function hermite_tracking(
    H::Union{Matrix, Vector}, 
    t::Number, 
    x::Vector, 
    r::Number, 
    A::AcbMatrix, 
    h::Number, 
    n::Int, 
    xprev::Vector, 
    vprev::Vector, 
    hprev::Number, 
    refinement_threshold::Number,
)
    Ft, x, r, A, v, h, radii = refine_step(H, t, x, r, A, h; threshold = refinement_threshold)
    X  = hermite_predictor(H, x, xprev, v, vprev, hprev)
    tm = krawczyk_operator_taylor_model(H, X, t, A, r)

    T = CCi("$radii +/- $radii")
    input = zeros(CCi, n + 1)
    input[n + 1] = T

    K = evaluate_matrix(tm, input)
    h = proceeding_step(h, CCi, n, tm, K)

    return x, v, h, X, r, A
end