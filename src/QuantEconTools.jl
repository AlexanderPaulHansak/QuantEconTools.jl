module QuantEconTools

export 
    # Markov functionality
    tauchen, rouwenhorst, stationary_distr,

    # Utility function
    AbstractUtility, LogUtility, CRRAUtility, derivative,

include("MarkovApprox.jl")
include("utility.jl")

end # module
