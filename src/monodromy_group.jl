export Vertex, vertex,
    Edge, edge,
    track_edge,
    search_point,
    track_complete_graph,
    complete_correspondences,
    neighbor,
    membership_test,
    p_compose,
    get_permutations,
    str_convert


struct Vertex
    base_point::Union{Vector{AcbFieldElem},Matrix{AcbFieldElem}}
    partial_sols::Vector{Union{Vector{AcbFieldElem},Matrix{AcbFieldElem}}}
    Edges::Vector{Any}
end

function Base.show(io::IO, x::Vertex)
    print(io, "Vertex")
end

struct Edge
    node1::Vertex
    node2::Vertex
    correspondence12::Vector{Tuple{Int64,Int64}}
    correspondence21::Vector{Tuple{Int64,Int64}}
end

function Base.show(io::IO, x::Edge)
    print(io, "Edge")
end

function vertex(p::Vector)
    Vertex(p, Matrix{AcbFieldElem}[], Vector{Edge}[])
end

function vertex(p::Vector, x::Vector)
    Vertex(p, x, Vector{Edge}[])
end

function edge(p::Vertex, q::Vertex)
    Edge(p, q, Tuple{Int64,Int64}[], Tuple{Int64,Int64}[])
end

function edge(
    p::Vertex,
    q::Vertex,
    c12::Vector{Tuple{Int64,Int64}},
    c21::Vector{Tuple{Int64,Int64}}
)
    Edge(p, q, c12, c21)
end

function track_edge(
    H::Union{Matrix, Vector}, 
    e::Edge, 
    from1to2::Bool, 
    r::Number;
    predictor = true
)
    if from1to2
        va, vb = e.node1, e.node2
        c12, c21 = e.correspondence12, e.correspondence21
    else
        va, vb = e.node2, e.node1
        c12, c21 = e.correspondence21, e.correspondence12
    end

    qa, qb = va.partial_sols, vb.partial_sols
    a, b = va.base_point, vb.base_point

    untracked_idx = setdiff(1:length(qa), map(t -> t[1], c12))
    println("→ tracking edge with $(length(untracked_idx)) paths")

    Fab = specified_system(a, b, H)
    n_iter = 1

    for i in qa[untracked_idx]
        x_idx = search_point(i, qa)
        x_idx === false && error("unknown candidate entered")

        println("\n🔁 tracking path #$(n_iter)")
        if predictor == false
            y = tracking_without_predictor(Fab, i, r)
        else
            y = track(Fab, i, r; show_display = true, refinement_threshold = 1/8)
        end
        n_iter += 1

        y_idx = search_point(y, qb)
        if y_idx === false
            push!(qb, y)
            y_idx = length(qb)
        end

        push!(c12, (x_idx, y_idx))
        push!(c21, (y_idx, x_idx))
    end

    va = vertex(a, qa)
    vb = vertex(b, qb)

    return from1to2 ? edge(va, vb, sort(c12), sort(c21)) :
                      edge(vb, va, sort(c21), sort(c12))
end

function search_point(res::Vector{AcbFieldElem}, p_list::Vector)
    n = length(p_list)
    best_idx, best_val = 0, Inf

    for i in 1:n 
        dist = maximum(map(abs ∘ max_int_norm, res - p_list[i]))
        if dist < best_val
            best_val = dist
            best_idx = i
        end
    end

    return best_val > 5e-3 ? false : Int64(best_idx)
end

function track_complete_graph(
    H::Union{Matrix,Vector},
    r::Number,
    vertices::Vector{Vertex},
    max_root_count::Int;
    predictor = true
)

    base_points = map(i -> i.base_point, vertices)
    rc = max_root_count
    edgs = Edge[]
    for i in 1:length(vertices)-1
        for j in i+1:length(vertices)
            e = edge(vertices[i], vertices[j])
            edgs = push!(edgs, e)

            node1_edges = vertices[i].Edges
            node2_edges = vertices[j].Edges
            node1_edges = unique(push!(node1_edges, e))
            node2_edges = unique(push!(node2_edges, e))

        end
    end
    npoints = map(i -> length(i.partial_sols), vertices)
    n_correspondences = map(i -> length(i.correspondence12), edgs)
    prev_n_correspondences = n_correspondences


    iter = 0
    while !all(==(max_root_count), n_correspondences)
        for i in edgs

            if sum(prev_n_correspondences) == sum(n_correspondences)
                iter = iter + 1
            else
                iter = 0
            end
            if iter > 10
                println("After 10 iterations, no known solution was found. Tracking is stopped.")
                return edgs
            end
            prev_n_correspondences = copy(n_correspondences)


            node1_idx = search_point(i.node1.base_point, base_points)
            node2_idx = search_point(i.node2.base_point, base_points)
            node1_sols = npoints[node1_idx]
            node2_sols = npoints[node2_idx]

            println("\r-----------------------------------------------------")
            println("\rstart node: $node1_idx : $node1_sols known solutions")
            println("\rtarg. node: $node2_idx : $node2_sols known solutions")
            if predictor == false
                i = track_edge(H, i, true, r; predictor=false)
            else
                i = track_edge(H, i, true, r)
            end

            println("\r-----------------------------------------------------")
            println("\rstart node: $node2_idx : $node2_sols known solutions")
            println("\rtarg. node: $node1_idx : $node1_sols known solutions")
            if predictor == false
                i = track_edge(H, i, false, r; predictor=false)
            else
                i = track_edge(H, i, false, r)
            end

            npoints = map(j -> length(j.partial_sols), vertices)
            n_correspondences = map(k -> length(k.correspondence12), edgs)

            println("\r# Correspondences per edge: $n_correspondences")
        end
    end

    edgs
end


function parameter_points(
    v1::Vertex,
    sz_p::Int,
    n_vertices::Int
)
    vertices = Vertex[v1]
    for i in 1:n_vertices-1
        v = AcbFieldElem[]
        for j in 1:sz_p
            r_unit_circle = exp(rand(Int8) * im)
            push!(v, CCi(real(r_unit_circle), imag(r_unit_circle)))
        end
        vertices = push!(vertices, vertex(v))
    end
    vertices
end

function complete_correspondences(
    rc::Int,
    E::Vector{Edge}
)
    reverse(deleteat!(E, findall(e ->
            length(unique(map(j -> j[2], e.correspondence12))) != rc || length(map(j -> j[2], e.correspondence21)) != rc, E)))
end

function neighbor(v::Vertex, e::Edge)
    if v == e.node1
        return e.node2
    elseif v == e.node2
        return e.node1
    else
        error("Edge is not incident at the given vertex.")
    end
end

function membership_test(
    v::Vertex,
    e::Edge,
    v_list::Vector{Vertex},
    e_list::Vector{Edge}
)
    u = neighbor(v, e)
    (u in v_list) == false && (e in e_list)
end

function p_compose(H1::Union{Vector{Tuple{Int64,Int64}},Vector{Pair{Int64,Int64}}}, H2::Vector{Pair{Int,Int}})

    l = length(H2)
    map(k -> k => sort(H1)[H2[k][2]][2], 1:l)
end

function get_permutations(
    rc::Number,
    E::Vector{Edge}
)

    id_perm = map(i -> i => i, 1:rc)
    EG = complete_correspondences(rc, E)
    VG = unique(vcat(map(i -> i.node1, EG), map(i -> i.node2, EG)))

    uncovered_v = deleteat!(VG, findall(i -> i == VG[1], VG))
    uncovered_e = EG

    T = []

    while length(uncovered_v) > 0
        v_list = [filter(v -> any(e -> membership_test(v, e, uncovered_v, uncovered_e), v.Edges), uncovered_v)[1]]
        if length(v_list) == 1
            v = v_list[1]

            e_list = [reverse(filter(e -> membership_test(v, e, uncovered_v, uncovered_e), v.Edges))[1]]
            if length(e_list) == 1
                e = e_list[1]
                T = push!(T, [v, e])
                uncovered_e = deleteat!(uncovered_e, findall(i -> i == e, uncovered_e))
            end
        end
        uncovered_v = deleteat!(uncovered_v, findall(i -> i == v, uncovered_v))
    end

    perms = Vector{Int64}[]
    for e in uncovered_e
        u = e.node1
        v = e.node2
        u_path = id_perm
        ind = findall(i -> u in i, T)
        while length(ind) > 0
            eu = T[ind[1]][2]
            if u == eu.node1
                u_path = p_compose(eu.correspondence12, u_path)
            else
                u_path = p_compose(sort(eu.correspondence21), u_path)
            end
            u = neighbor(u, eu)
            ind = findall(i -> u in i, T)
        end
        v_path = id_perm
        ind = findall(i -> v in i, T)
        while length(ind) > 0
            ev = T[ind[1]][2]
            if v == ev.node1
                v_path = p_compose(ev.correspondence12, v_path)
            else
                v_path = p_compose(sort(ev.correspondence21), v_path)
            end
            v = neighbor(v, ev)
            ind = findall(i -> v in i, T)
        end
        perms = push!(perms, map(i -> i[2], p_compose(v_path, p_compose(e.correspondence12, sort(map(i -> reverse(i), u_path))))))
    end

    perms
end

function safe_path(file_name::String)
    if startswith(file_name, "~")
        return expanduser(file_name)
    elseif isabspath(file_name)
        return file_name
    else
        return joinpath(pwd(), file_name)
    end
end

function str_convert(
    l::Vector{Vector{Int64}},
    file_name::String,
    group_name::String
)
    path = safe_path("$file_name.txt")
    open(path, "w") do file ## change the directory properly
        iter = 0
        for i in l
            write(file, "p$iter:= PermList($i);\n")
            iter = iter + 1
        end

        write(file, "$group_name:=Group(")
        for i in 0:length(l)-1
            if i < length(l) - 1
                write(file, "p$i,")
            else
                write(file, "p$i")
            end
        end
        write(file, ");")
    end
end



@gap("GaloisWidth := function(G)
  local X, M, C, phi;
  if IsTrivial(G) then return 1;
  elif IsNaturalSymmetricGroup(G) or IsNaturalAlternatingGroup(G) then
    X := OrbitsDomain(G)[1];
    if Length(X) = 4 then return 3;
    else return Length(X);
    fi;
  elif IsCyclic(G) then return Maximum(Factors(Order(G)));
  elif not IsTransitive(G) then 
    return Maximum(List(Orbits(G), 
      O -> GaloisWidth(Image(ActionHomomorphism(G,O)))
    ));
  else
    X := OrbitsDomain(G)[1];
    if not IsPrimitive(G) then
      phi := ActionHomomorphism(G,Blocks(G, X),OnSets);
      return Maximum(GaloisWidth(Kernel(phi)), GaloisWidth(Image(phi)));
    elif IsSimple(G) then
      M := List(ConjugacyClassesMaximalSubgroups(G), H -> Representative(H));
      return Minimum(List(M, H -> Order(G)/Order(H)));
    else
      C := CompositionSeries(G);
      return Maximum(List([1..Length(C)-1], 
        i -> GaloisWidth(C[i]/C[i+1])
      ));
    fi;
  fi;
end;
")