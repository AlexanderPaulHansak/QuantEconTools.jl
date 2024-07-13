
using Test, QuantEconTools # This load both the test suite and our Package


# Test some markov functionalities
@test all(abs.(inv_distr(tauchen(0.9,0.2,3,3).Π).InvD' - inv_distr(tauchen(0.9,0.2,3,3).Π).InvD'*tauchen(0.9,0.2,3,3).Π) .<= 1e-10)
@test all(abs.(inv_distr(rouwenhorst(0.9,0.2,3).Π).InvD' - inv_distr(rouwenhorst(0.9,0.2,3).Π).InvD'*rouwenhorst(0.9,0.2,3).Π) .<= 1e-10)