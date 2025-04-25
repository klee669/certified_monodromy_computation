include("../../../../src/certified_monodromy_computation.jl")


@setupfield begin
    AcbField(1024)
    (x,)
    (η,)
    (t,)
    (c0,c1,c2,c3,c4,c5,c6,c7,c8)
end

CCi = _CCi
r = .1;

F = hcat([c8*x^8+c7*x^7+c6*x^6+c5*x^5+c4*x^4+c3*x^3+c2*x^2+c1*x+c0])
bp = [CCi(1/2), CCi(9/7), CCi(9/7), CCi(45/7), CCi(45/56), CCi(-2/3), CCi(1/3), CCi(3/4), CCi(3/7)]
x = [CCi(.0720838,-.481488)]
v1 = vertex(bp,[x])
n_nodes = 7;

@gap("Read(\"~/Documents/GitHub/certified_monodromy_comp/examples/univariate_examples/sum_of_even_degrees/even_degree_sums_8.txt\");")
@gap("G;")
@gap("StructureDescription(G);") # ((((C2 x C2 x C2) : (C2 x C2)) : C3) : C2) : C2
@gap("GaloisWidth(G);") #3


path = safe_path("~/Documents/GitHub/certified_monodromy_comp/examples/univariate_examples/sum_of_even_degrees/results_$(n_nodes)_nodes.txt")
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
            vs = parameter_points(v1, 9, n_nodes)
            edges = track_complete_graph(F, r, vs,8)
            if length(edges[1].correspondence12) != 8
                fail_correspondence_count = fail_correspondence_count+1;
            end

            perms=get_permutations(length(edges[1].correspondence12),edges)
            str_convert(perms, "~/Documents/GitHub/certified_monodromy_comp/examples/univariate_examples/sum_of_even_degrees/even_degree_sums_8_$(n_nodes)nodes", "H")

            filename = expanduser("~/Documents/GitHub/certified_monodromy_comp/examples/univariate_examples/sum_of_even_degrees/even_degree_sums_8_$(n_nodes)nodes.txt")
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
                write(file, "⚠️ Error at i=$i: $(e)\n\n")
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


