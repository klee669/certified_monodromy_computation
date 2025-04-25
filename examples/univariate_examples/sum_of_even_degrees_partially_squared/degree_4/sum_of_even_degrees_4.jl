include("../../../../src/certified_monodromy_computation.jl")


@setupfield begin
    AcbField()
    (x,)
    (η,)
    (t,)
    (c1,c2,c3)
end

CCi = _CCi
r = .1;

F = hcat([c3^2*x^4+c2*x^2+c1^2])
bp = [CCi(10/3), CCi(1/3), CCi(5/4)]
x = [CCi(.933811,.784519)]
v1 = vertex(bp,[x])
n_nodes = 7;

@gap("Read(\"~/Documents/GitHub/certified_monodromy_comp/examples/univariate_examples/sum_of_even_degrees_partially_squared/degree_4/even_degree_sums_4.txt\");")
@gap("G;")
@gap("StructureDescription(G);") # C2 x C2
@gap("GaloisWidth(G);") #2

for n_nodes in 3:6
    path = safe_path("~/Documents/GitHub/certified_monodromy_comp/examples/univariate_examples/sum_of_even_degrees_partially_squared/degree_4/results_degree_4_$(n_nodes)_nodes.txt")

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
            vs = parameter_points(v1, 4, n_nodes)
            edges = track_complete_graph(F, r, vs,6)
            if length(edges[1].correspondence12) != 6
                fail_correspondence_count = fail_correspondence_count+1;
            end

            perms=get_permutations(length(edges[1].correspondence12),edges)
            str_convert(perms, "~/Documents/GitHub/certified_monodromy_comp/examples/univariate_examples/sum_of_even_degrees_partially_squared/degree_4/even_degree_sums_4_$(n_nodes)nodes", "H")
    
            filename = expanduser("~/Documents/GitHub/certified_monodromy_comp/examples/univariate_examples/sum_of_even_degrees_partially_squared/degree_4/even_degree_sums_4_$(n_nodes)nodes.txt")
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
            #                write(file, "⚠️ Error at i=$i: $(e)\n\n")
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


end