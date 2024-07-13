
using Test, QuantEconTools # This load both the test suite and our Package

out = plusTwo(3)

@test out == 5  