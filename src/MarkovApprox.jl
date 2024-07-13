#=
Various routines to discretize AR(1) processes

@author : Alexander Hansak <alexander.hansak@cerge-ei.cz>

@date : 2024-04-10 23:55:05

=#


import Distributions: cdf, Normal
import LinearAlgebra: eigvecs, eigvals


"""
    MarkovApprox(ρ,σ,m,N,Method)

Returns transition matrix, invariant distribution and discretized vector for approximating AR(1) process with either Tauchen (1986) or Rouwenhorst (1995)

# Examples
```julia-repl
julia> MarkovApprox(0.5,1,3,3,0).s
-3.4641016151377544:3.4641016151377544:3.4641016151377544

julia> MarkovApprox(0.5,1,3,3,0).Π
3×3 Matrix{Float64}:
 0.5          0.499734  0.000266003
 0.0416323    0.916735  0.0416323
 0.000266003  0.499734  0.5

julia> MarkovApprox(0.5,1,3,3,0).InvD
3-element Vector{Float64}:
 0.3082000863921907
 0.3835998272156185
 0.3082000863921907
```
"""
function tauchen(ρ,σ,m,N)
    #= 
    Implements Tauchen's (1986) or Rouwenhorst's (1995)
    method to discretize the first-order 
    autoregressive process 
                               y_t = ρ * y_(t-1) + u_t
    with a Markov chain.
    
    u_t is Gaussian white noise with standard deviation σ.
    
    ----------------------------------------------------------------------------------------
    INPUT:  ρ   autocorrelations coefficient
            σ   tandard deviation of Gaussian white noise
            m   width of discretized state space (ymax=m*vary ymin=-m*vary, Tauchen uses m=3) 
            N   number of possible states to approximate y_t (usually N=9 should be fine)
    
    OUTPUT: Π   the transition matrix of the Markov chain
            s   the discretized state space of y_t
    -----------------------------------------------------------------------------------------
    =#

    # Check if abs(ρ)<1 
    if abs(ρ)>=1
        error("The persistence parameter, ρ, must be less than one in absolute value.")
    end
    
    Π = zeros(N,N);
    dist = Normal(0,σ);   # we use normal distribution with mu=0 for u_t
    
    
    # discretize the state space
    stvy = sqrt(σ^2/(1-ρ^2));   # standard deviation of y_t
    ymax = m*stvy;                    # upper boundary of state space
    ymin = -ymax;                     # lower boundary of state space
    w = (ymax-ymin)/(N-1);            # length of interval 
    s = ymin:w:ymax;                  # the discretized state space
    
    # calculate the transition matrix
    for j=1:N
        for k=2:N-1
            Π[j,k]= cdf(dist,s[k]-ρ*s[j]+w/2) - cdf(dist,s[k]-ρ*s[j]-w/2);
        end
        Π[j,1] = cdf(dist,s[1]-ρ*s[j]+w/2);
        Π[j,N] = 1 - cdf(dist,s[N]-ρ*s[j]-w/2);
    end

    return (; Π,s)

    # Check for correct Transition matrix
    if any(abs.(sum(Π,dims=2).-1) .> 1e-10)
        str = findall(>(1e-10),vec(abs.(sum(Π,dims=2).-1)));  # find rows not adding up to one
        println("error in transition matrix")
        println("rows $str do not sum to one")
    end
    
end
    
    

function rouwenhorst(ρ,σ,N)
# Use Rouwenhorst's (1995) method

    # discretize the state space
    stvy = sqrt(σ^2/(1-ρ^2)); # standard deviation of y_t
    ymax = sqrt(N-1)*stvy;             # upper boundary of state space
    ymin = -ymax;                      # lower boundary of state space
    w = (ymax-ymin)/(N-1);             # length of interval 
    s = ymin:w:ymax;                   # the discretized state space
    
    # calculate the transition matrix
    p       = (1+ρ)/2;
    q       = p;
    Omega0  = [p 1-p; 1-q q];
    
    for i = 3:N
        zero        = zeros(i-1,1);
        OmegaNew    = p.*[Omega0 zero; zero' 0] + (1-p).*[zero Omega0; 0 zero']+(1-q).*[zero' 0; Omega0 zero] + q.*[0 zero'; zero Omega0];
        Omega0      = OmegaNew;
        Omega0[2:end-1,:] = Omega0[2:end-1,:]./2;
    end
    
    Π = Omega0;
    
    # Check for correct transition matrix
    if any(abs.(sum(Π,dims=2).-1) .> 1e-10)
        str = findall(>(1e-10),vec(abs.(sum(Π,dims=2).-1)));  # find rows not adding up to one
        println("error in transition matrix")
        println("rows $str do not sum to one")
    end

return (; Π,s)

end

function inv_distr(Π)
    
    # Calculate invariant distribution over states
    V = eigvals(Π');
    Σ = eigvecs(Π');
    InvD = [];
    for i in eachindex(V)
        if (abs(V[i]-1)<1e-14)
            InvD = Σ[:,i]/sum(Σ[:,i]);
        end
    end
    return (; InvD)
end

function inv_distr2(Π)
    
    # Calculate invariant distribution over states
    InvD = Π^100

    return(; InvD)
end
       
       
    
