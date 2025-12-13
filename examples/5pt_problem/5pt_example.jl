include("../../src/certified_monodromy_computation.jl")

# polynomial rings
@setupfield begin
    AcbField()
    (X,Y,Z)
    (η,)
    (t,)
    (x11,x21,x31,x41,x51,x12,x22,x32,x42,x52,y11,y21,y31,y41,y51,y12,y22,y32,y42,y52)
end
CCi = _CCi

include("eqs.jl")
include("pts.jl")

#base_point
bp = [CCi(.831187,-.555993), CCi(-.077487,.996993), CCi(-.031994,-.999488), CCi(.019587,.999808), CCi(.949892,-.312579), CCi(-.881637,.471928), CCi(.203765,.97902), CCi(-.819274,.573402), CCi(-.973777,-.227505), CCi(-.972427,-.233208), CCi(.192661,-.981265), CCi(-.856378,-.516349), CCi(.296842,-.954927), CCi(.525122,-.851027), CCi(-.999761,-.021852), CCi(-.679705,-.733486), CCi(-.583596,.812044), CCi(.875249,.483673), CCi(.43059,-.902548), CCi(.438152,.898901)]
x0 = p_list[1]

v1 = vertex(bp,[x0])
vs = parameter_points(v1, 20, 4)
r = .1;
edges = track_complete_graph(F, r, vs,20)

perms=get_permutations(20,edges)
str_convert(perms, "5pt_perm_list_20")

using GAP
@gap("Read(\"~/Documents/GitHub/certified_monodromy_comp/examples/5pt_problem/5pt_perm_list_20.txt\");")
@gap("G;")
@gap("StructureDescription(G);") # (C2 x C2 x C2 x C2 x C2 x C2 x C2 x C2 x C2) : S10
@gap("GaloisWidth(G);") #10


