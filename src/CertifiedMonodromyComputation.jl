module CertifiedMonodromyComputation

using Nemo
using AbstractAlgebra
using LinearAlgebra
using MultivariatePolynomials
using GAP

export @setupfield

export straight_line_homotopy, specified_system, track

export track_complete_graph, get_permutations, str_convert

# Source Code Include

# [Internals] -- helpers and utilities
include("internals/interval_arithmetic.jl")
include("internals/linear_algebra.jl")
include("internals/moore_box.jl")       
include("internals/krawczyk.jl")        
include("internals/predictors.jl")
include("internals/tracking_modules.jl") 

# [Core] -- main functionalities
include("poly_setup.jl")
include("homotopy.jl")    
include("tracking.jl")    
include("monodromy.jl")   

end # module