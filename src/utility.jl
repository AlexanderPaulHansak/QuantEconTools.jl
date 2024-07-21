abstract type AbstractUtility end

#
# Separable utility
#

# Consumption utility

@doc doc"""
Type used to evaluate log utility. Log utility takes the form

u(c) = \log(c)

Additionally, this code assumes that if c < 1e-10 then

u(c) = log(1e-10) + 1e10*(c - 1e-10)

"""
struct LogUtility <: AbstractUtility
    ξ::Float64
end

LogUtility() = LogUtility(1.0)

(u::LogUtility)(c::Float64) =
    c > 1e-10 ? u.ξ*log(c) : u.ξ*(log(1e-10) + 1e10*(c - 1e-10))
derivative(u::LogUtility, c::Float64) =
    c > 1e-10 ? u.ξ / c : u.ξ*1e10

"""
Type used to evaluate CRRA utility. CRRA utility takes the form

u(c) = c^(1 - γ) / (1 - γ)

Additionally, this code assumes that if c < 1e-10 then

u(c) = (1e-10^(1 - γ) / (1 - γ) + 1e-10^(-γ) * (c - 1e-10))
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

(u::CRRAUtility)(c::Float64) =
    c > 1e-10 ?
           (c^(1.0 - u.γ) - 1.0) / (1.0 - u.γ) :
           ((1e-10^(1.0 - u.γ) - 1.0) / (1.0 - u.γ) + 1e-10^(-u.γ)*(c - 1e-10))
derivative(u::CRRAUtility, c::Float64) =
    c > 1e-10 ? c^(-u.γ) : 1e-10^(-u.γ)
