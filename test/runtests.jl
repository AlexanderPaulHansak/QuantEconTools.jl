
using Test, QuantEconTools # This load both the test suite and our Package

out = tauchen(0.9,0.2,3,3).s

@test out == 5  