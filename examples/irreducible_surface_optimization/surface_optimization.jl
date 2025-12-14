using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
#Pkg.instantiate()

using Nemo
using AbstractAlgebra
using CertifiedMonodromyComputation
using GAP

# ------------------------------------------------------------------------------
# 1. Setup Polynomial Ring & System
# ------------------------------------------------------------------------------
# @setupfield to define variables and fields
@monodromy_setup begin
    vars = (x,y,z,λ)
    params = (u1,u2,u3)
end
const CCi = _CCi

F = [2*(x-u1)-6*λ*(x^2+y^2)^2*x 2*(y-u2)-6*λ*(x^2+y^2)^2*y 2*(z-u3)+4*λ*z^3 0*u1+z^4-(x^2+y^2)^3]
vars = [x y z λ]
pars = [u1 u2 u3]

# ------------------------------------------------------------------------------
# 2. Setup Vertices & Initial Points
# ------------------------------------------------------------------------------
# define base_point 
bp = [CCi(.09868,.675389), CCi(.423238,.713082), CCi(.592351,.144969)]
x0 = [CCi(1.23836,-.422501), CCi(1.19574,-1.0474), CCi(2.08916,1.85256), CCi(-.0126777,.0505892)]

v1 = vertex(bp,[x0])
vs = parameter_points(v1, 3, 4)

# ------------------------------------------------------------------------------
# 3. Solve Monodromy (Tracking)
# ------------------------------------------------------------------------------
println("Starting Monodromy Tracking...")

edges = solve_monodromy(F, vs; max_roots=8)

# ------------------------------------------------------------------------------
# 4. GAP Group Construction
# ------------------------------------------------------------------------------
println("Building GAP Group...")

G = build_gap_group(8, edges)


# ------------------------------------------------------------------------------
# 5. Analysis (GAP)
# ------------------------------------------------------------------------------
if G !== nothing
    println("Structure Description:")
    println(GAP.Globals.StructureDescription(G)) # (S4 x S4) : C2
    println("Galois Width:")
    gw = galois_width(G)
    println(gw) # 3
end

