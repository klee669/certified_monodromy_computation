export track, tracking_without_predictor

# tracking without predictor
function tracking_without_predictor(H, x; r = .1, iterations_count = false)
    try
        ring = parent(H[1]);
        coeff_ring = base_ring(H[1]);
        CCi = coefficient_ring(coeff_ring)
        t = 0;
        h = 1;
        n = length(x);
        G = evaluate_matrix(Matrix(transpose(hcat(H))), 0)
        A = jacobian_inverse(G, x)

        iter = 0;
        while t < 1
            yield()
            
            rt = round(t, digits = 10)
            print("\r(h, progress t): ($h,$rt)")

            Ft = evaluate_matrix(Matrix(transpose(hcat(H))), t);
            x,r,A = refine_moore_box(Ft, x, r, A, 1/8);
            h = 2*h;
            radii = h/2;


            midt = t+h/2;
            T = CCi("$midt +/- $radii");
            FT = evaluate_matrix(Matrix(transpose(hcat(H))), t+h);
            while krawczyk_test(FT, x, r, A, 7/8) == false
                h = 1/2 * h;
                midt = t+h/2;
                radii = h/2;
        
                T = CCi("$midt +/- $radii");
                FT = evaluate_matrix(Matrix(transpose(hcat(H))), t+h);
            end
            t = max_int_norm(T);
            iter = iter+1;
        end

        Ft = evaluate_matrix(Matrix(transpose(hcat(H))), 1);
        x,r,A = refine_moore_box(Ft, x, r, A, 1/8);
        
        if iterations_count
            return x, iter
        else
            return x
        end

    catch e
        if isa(e, InterruptException)
            println("\n[Warning] Tracking interrupted manually.")
            return nothing
        else
            rethrow(e)
        end
    end
end


"""
    track(H, x, r; options...)

Track a solution path. Returns `nothing` if interrupted.
"""
function track(
    H::Union{Matrix, Vector},
    x::Vector{AcbFieldElem};
    r = 0.1, # initial radius
    show_display = true,
    refinement_threshold = 1/8,
    predictor = "hermitian",
    iterations_count = false,
    tracking = "truncate",
    projective = false,
)
    try
        if predictor == "without_predictor"
            return tracking_without_predictor(H, x)
        end

        ## --- Initialization --------------------------------------------------
        coeff_ring  = base_ring(H[1])
        CCi         = coefficient_ring(coeff_ring)
        t           = 0
        h           = 1 / 10
        iter        = 0
        tprev       = 0
        n           = length(x)
        max_deg_H   = max_degree(hcat(H))

        G = evaluate_matrix(Matrix(transpose(hcat(H))), 0)
        A = jacobian_inverse(G, x)

        ## --- First Step (Linear predictor) ----------------------------------
        x, v, h, X, r, A = linear_tracking(H, t, x, r, A, h, n, refinement_threshold)

        xprev  = x
        vprev  = v
        hprev  = h
        t     += h

        ## --- Main Loop -------------------------------------------------------
        while t < 1
            yield()
            rt = round(t, digits = 10)

            H, x, v, h, X, r, A = hermite_tracking(
                H, max_deg_H, t, x, r, A, h, n, xprev, vprev, hprev, refinement_threshold;
                tracking, projective
            )

            xprev  = x
            vprev  = v
            hprev  = h
            t     += h

            input        = zeros(CCi, n + 1)
            input[n + 1] = CCi(h)

            x = midpoint_complex_box(Matrix(evaluate_matrix(matrix(X), input))[1:n])
            iter += 1

            show_display && print("\rprogress t: $rt")
        end

        show_display && print("\r ")

        ## --- Final Refinement ------------------------------------------------
        Ft, x, r, A, v, h, radii = refine_step(H, 1, x, r, A, h; threshold = 1 / 100)
        
        if iterations_count
            return x, iter
        else
            return x
        end

    catch e
        if isa(e, InterruptException)
            println("\n[Warning] Tracking interrupted manually.")
            return nothing 
        else
            rethrow(e)
        end
    end
end