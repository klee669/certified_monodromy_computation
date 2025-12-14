export speed_vector,
       linear_predictor,
       hermite_predictor


# --------------------------------------------------------------------------
# Compute the speed vector for the linear predictor
function speed_vector(
    H::Union{Matrix, Vector}, 
    x::Vector{AcbFieldElem}, 
    tval::Number, 
    A::AcbMatrix,
)
    ring = base_ring(H[1])
    n    = length(H)

    # Evaluate dH/dt symbolically at t
    result = evaluate_matrix(Matrix(transpose(hcat(H))),tval);
    for i = 1:n
        result[i]=Nemo.evaluate(derivative(H[i]),tval);
    end
    result = evaluate_matrix(result, x);

    midpoint_complex_box(-A*transpose(hcat(result)))
end


# --------------------------------------------------------------------------
# Linear predictor step: x + η·v
function linear_predictor(
    H::Union{Matrix, Vector}, 
    v::Vector{AcbFieldElem}, 
    x::Vector{AcbFieldElem},
)
    eR    = base_ring(H[1])
    η     = gens(eR)[end]
    delta = [vi * η for vi in v]

    return x + delta
end


# --------------------------------------------------------------------------
# Hermite predictor (uses previous step data)
function hermite_predictor(
    H::Union{Matrix, Vector}, 
    x::Vector, 
    xprev::Vector, 
    v::Vector, 
    vprev::Vector, 
    hprev::Number,
)
    eR = base_ring(H[1])
    η  = gens(eR)[end]
    n  = length(v)

    result = Vector{typeof(η)}(undef, n)

    for i in 1:n
        dη   = v[i] * η
        η²   = η^2
        η³   = η^3

        dx   = x[i] - xprev[i]
        dv   = v[i] - vprev[i]
        h²   = hprev^2
        h³   = hprev^3

        t2 = (3 * v[i] / hprev - dv / hprev - 3 * dx / h²) * η²
        t3 = (2 * v[i] / h²   - dv / h²   - 2 * dx / h³) * η³

        result[i] = dη + t2 + t3
    end

    return x + result
end
