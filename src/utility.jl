abstract type AbstractUtility end

#
# Separable utility
#

# Consumption utility

"""
Type used to evaluate log utility. Log utility takes the form

u(c) = log(c)

Additionally, this code assumes that if c < 1e-10 then

u(c) = -10e18

"""
struct LogUtility <: AbstractUtility end

(u::LogUtility)(c::Float64) =
    c > 1e-10 ? log(c) : -1e18
derivative(u::LogUtility, c::Float64) =
    c > 1e-10 ? 1 / c : 1e10

"""
Type used to evaluate CRRA utility. CRRA utility takes the form

u(c) = c^(1 - σ) / (1 - σ)

Additionally, this code assumes that if c < 0 then

u(c) = -1e18
"""
struct CRRAUtility <: AbstractUtility
    γ::Float64

    function CRRAUtility(γ)
        if abs(γ - 1.0) < 1e-8
            error("Your value for γ is very close to 1... Consider using LogUtility")
        end

        return new(γ)
    end
end

(u::CRRAUtility)(c::Real) =
    c > 0 ?
           c^(1.0 - u.γ) / (1.0 - u.γ) : -1e18
derivative(u::CRRAUtility, c::Float64) =
    c > 1e-10 ? c^(-u.γ) : 1e-10^(-u.γ)


# Labor Utility

"""
Type used to evaluate constant Frisch elasticity (CFE) utility. CFE
utility takes the form

v(h) =  χ h^(1 + 1/γ) / (1 + 1/γ)

Code assumes input h ∈ [0,1], will not correct for different values.
"""
struct CFEUtility <: AbstractUtility
    γ::Float64
    χ::Float64

    function CFEUtility(γ, χ=1.0)
        return new(γ, χ)
    end
end

(u::CFEUtility)(h::Float64) =
    -u.χ * h^(1.0 + 1.0/u.γ)/(1.0 + 1.0/u.γ) 
derivative(u::CFEUtility, h::Float64) =
    -u.χ * h^(1.0/u.γ) 

