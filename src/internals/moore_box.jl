export refine_moore_box

# refining the Moore box
function refine_moore_box(
    system::Union{Matrix,AbstractAlgebra.Generic.MatSpaceElem}, 
    point::Vector{AcbFieldElem}, 
    r::Number,
    A::AcbMatrix,
    ρ::Number
)
    y = point;
    CR = parent(point[1]);
    n = size(A)[1];
    while krawczyk_test(system, y, r, A, ρ) == false
        d = A * transpose(evaluate_matrix(system, y));
        if max_norm(d) <= (1/64)*ρ*r
            r = (1/2)*r;
        else
            y = midpoint_complex_box(y-d[:,1]);
        end
        A = jacobian_inverse(system, y);
    end
    while 2*r <= 1 && krawczyk_test(system, point, 2*r, A, ρ)
        r = 2*r;
    end

    [y, r, A]
end