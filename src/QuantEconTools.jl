module QuantEconTools

export 
    # Markov functionality
    tauchen, rouwenhorst, stationary_distr

include("MarkovApprox.jl")
include("utility.jl")

end # end module
