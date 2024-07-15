#=
Various routines for discretizing and working with AR(1) processes

@author : Alexander Hansak <alexander.hansak@cerge-ei.cz>

@date : 2024-07-10 14:14:05

=#


import Distributions: cdf, Normal
import LinearAlgebra: I

@inline is_stochastic(P) = maximum(abs, sum(P, dims = 2) .- 1) < 5e-15 ? true : false

"""
Implements Tauchen's (1986) method for approximating AR(1) process with finite markov chain

The process follows

```math
    y_t = ρ y_{t-1} + ε_t
```

where ``ε_t ∼ N(0,σ^2)``

##### Arguments

- `ρ::Real` : Persistence parameter in AR(1) process
- `σ::Real` : Standard deviation of random component of AR(1) process
- `m::Real(3)` : The number of standard deviations the process should span to each side
- `N::Integer(9)`: Number of states in markov process

##### Returns

- `P::Matrix` : Transition matrix where P(i,j) denotes probability of going from state i to state j
- `s::StepRangeLen` : Discretized grid of AR(1) process according to Markov chain states
"""
function tauchen(ρ::Real,σ::Real,m::Real=3,N::Integer=9)
    
    # Check if abs(ρ)<1 
    abs(ρ)>=1 &&  
        throw(ArgumentError("Persistence parameter ρ must be smaller than 1"))
    
    P = zeros(N,N);
    norm_cdf(x) = cdf(Normal(0,σ),x)  # we use normal distribution with mean zero

    # discretize the state space
    stvy = sqrt(σ^2/(1-ρ^2))       # standard deviation of y_t
    ymax = m*stvy                  # upper boundary of state space
    ymin = -ymax                   # lower boundary of state space
    w = (ymax-ymin)/(N-1)          # length of interval 
    s = range(ymin,ymax,N)         # the discretized state space
    
    # calculate the transition matrix
    for row = 1:N

        for col = 2:N-1
            P[row,col]= norm_cdf(s[col]-ρ*s[row]+w/2) -
                        norm_cdf(s[col]-ρ*s[row]-w/2)
        end

        P[row,1] = norm_cdf(s[1]-ρ*s[row]+w/2);
        P[row,N] = 1 - norm_cdf(s[N]-ρ*s[row]-w/2)
    end

    return (; P,s)

    # Check for correct Transition matrix
    !is_stochastic(P) &&
        println("Problem in Transition matrix: Rows do not sum to 1!")
    
end
    
    
"""
Implements Rouwenhorst's (1995) method for approximating AR(1) process with finite markov chain

The process follows

```math
    y_t = ρ y_{t-1} + ε_t
```

where ``ε_t ∼ N (0,σ^2)``

##### Arguments

- `ρ::Real` : Persistence parameter in AR(1) process
- `σ::Real` : Standard deviation of random component of AR(1) process
- `N::Integer`: Number of states in markov process

##### Returns

- `P::Matrix` : Transition matrix where P(i,j) denotes probability of going from state i to state j
- `s::StepRangeLen` : Discretized grid of AR(1) process according to Markov chain states
"""
function rouwenhorst(ρ::Real,σ::Real,N::Integer)

    # Check if abs(ρ)<1 
    abs(ρ)>=1 &&  
        throw(ArgumentError("Persistence parameter ρ must be smaller than 1"))

    # discretize the state space
    stvy = σ/sqrt(1-ρ^2)       # standard deviation of y_t
    ymax = sqrt(N-1)*stvy      # upper boundary of state space
    s = range(-ymax, ymax,N)   # the discretized state space
    
    # The initial transition matrix for N=2
    p  = (1+ρ)/2
    q  = p
    P  = [p 1-p; 1-q q]
    
    # Recursively calculate the transition matrix for n > 2
    for n = 2:N-1
        P_old = P
        P = fill(0.,(n+1,n+1))
        P[1:n, 1:n]         = p*P_old
        P[1:n, 2:(n+1)]     .+= (1-p)*P_old
        P[2:(n+1), 1:n]     .+= (1-q)*P_old
        P[2:(n+1), 2:(n+1)] .+= q*P_old

        P[2:end-1,:] ./= 2
    end
    
    # Check for correct Transition matrix
    !is_stochastic(P) &&
        println("Problem in Transition matrix: Rows do not sum to 1!")

    return (; P,s)

end



"""
Computes the stationary distribution of an irreducible Markov
transition matrix P (stochastic matrix where rows sum to one).

Specifically, computes vector μ with 

```math
    μ * P = μ
```

and `` sum(μ) = 1 ``

Method: Comupute normalized eigenvector corresponding to left-eigenvector λ=1

##### Arguments

- `P::Matrix{T}` : Stochastic matrix of size n x n, elements must be non-negative and rows must sum to 1.

##### Returns

- `μ::Vector{T}` : stationary distribution of ``P``.
"""
function stationary_distr(P::Matrix{<:Real})

    # Check for correct Transition matrix
    !is_stochastic(P) && 
        throw(ArgumentError("Non-stochastic transition matrix provided"))
    
    # find an eigenvector x corresponding to eigenvalue λ = 1,
    # normalized so that the first component is 1
    A = P'
    μ = [1; (I - A[2:end,2:end]) \ Vector(A[2:end,1])]
    μ /= sum(μ)
    return μ
end


       
    
