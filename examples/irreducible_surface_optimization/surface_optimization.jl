include("../../src/certified_monodromy_computation.jl")

# polynomial rings
@setupfield begin
    AcbField()
    (x,y,z,λ)
    (η,)
    (t,)
    (u1,u2,u3)
end
CCi = _CCi

F = [2*(x-u1)-6*λ*(x^2+y^2)^2*x 2*(y-u2)-6*λ*(x^2+y^2)^2*y 2*(z-u3)+4*λ*z^3 0*u1+z^4-(x^2+y^2)^3]
vars = [x y z λ]
pars = [u1 u2 u3]

bp = [CCi(.09868,.675389), CCi(.423238,.713082), CCi(.592351,.144969)]
x0 = [CCi(1.23836,-.422501), CCi(1.19574,-1.0474), CCi(2.08916,1.85256), CCi(-.0126777,.0505892)]

v1 = vertex(bp,[x0])
vs = parameter_points(v1, 3, 4)
r = .1
edges = track_complete_graph(F, r, vs,8)

perms=get_permutations(8,edges)
str_convert(perms, "surface_optimization_8")

using GAP
@gap("Read(\"~/Documents/GitHub/certified_monodromy_comp/examples/irreducible_surface_optimization/surface_optimization_8.txt\");")
@gap("G;")
@gap("StructureDescription(G);") # (S4 x S4) : C2
@gap("GaloisWidth(G);") #3