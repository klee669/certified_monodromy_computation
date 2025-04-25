include("../../../../src/certified_monodromy_computation.jl")


@setupfield begin
    AcbField(1024)
    (x,)
    (η,)
    (t,)
    (c0,c1,c2,c3,c4)
end

CCi = _CCi
r = .1;

F = hcat([c4*x^4+c3*x^3+c2*x^2+c1*x+c0])
bp = [CCi(1/2), CCi(9/7), CCi(9/7), CCi(45/7), CCi(45/56)]
x = [CCi(.0720838,-.481488)]
v1 = vertex(bp,[x])
n_nodes = 6;

@gap("Read(\"~/Documents/GitHub/certified_monodromy_comp/examples/univariate_examples/random_polynomials/degree_4/random_polynomials_4.txt\");")
@gap("G;")
@gap("StructureDescription(G);") # S4
@gap("GaloisWidth(G);") # 4


path = safe_path("~/Documents/GitHub/certified_monodromy_comp/examples/univariate_examples/random_polynomials/degree_4/results_$(n_nodes)_nodes.txt")
using GAP
gw_counts = Dict{Int, Int}()

#@gap("LoadPackage(\"sonata\");")
open(path, "w") do file
    false_count = 0;
    fail_correspondence_count = 0;
    tracking_error_count = 0;
    for i in 1:100
        try
            v1 = vertex(bp,[x])
            vs = parameter_points(v1, 5, n_nodes)
            edges = track_complete_graph(F, r, vs,4)
            if length(edges[1].correspondence12) != 4
                fail_correspondence_count = fail_correspondence_count+1;
            end

            perms=get_permutations(length(edges[1].correspondence12),edges)
            str_convert(perms, "~/Documents/GitHub/certified_monodromy_comp/examples/univariate_examples/random_polynomials/degree_4/random_polynomials_4_$(n_nodes)nodes", "H")

            filename = expanduser("~/Documents/GitHub/certified_monodromy_comp/examples/univariate_examples/random_polynomials/degree_4/random_polynomials_4_$(n_nodes)nodes.txt")
            cmd = string("Read(\"", filename, "\");")
            GAP.evalstr(cmd)
            A=@gap("StructureDescription(H);") 
            println(A);
            tf=string(@gap("IsomorphismGroups(G,H);"))
            if tf == "fail"
                gw = @gap("GaloisWidth(H);")
                gw_counts[gw] = get(gw_counts, gw, 0) + 1
                false_count = false_count+1;
                tf = "false";
            else
                tf = "true";
            end
            write(file, "G#$i description: $A\n")
            write(file, "    isomorphic?: $tf\n")
            if tf == "false"
                write(file, "    Galois width: $gw\n")
            end
            write(file, "\n")
            flush(file)

            catch e
                println("⚠️ Error at i=$i: $(e)")
                #                    write(file, "⚠️ Error at i=$i: $(e)\n\n")
                i = i-1;
            continue  
        end
            
    end
    write(file, "number of incomplete correspondences: $fail_correspondence_count\n")
    write(file, "false count: $false_count\n")

    write(file, "Galois width counts:\n")
    for (k,v) in sort(collect(gw_counts))
        write(file, "gw = $k occurred $v times\n")
    end
end


