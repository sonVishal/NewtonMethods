#!/usr/bin/env julia

# Start Test Script
using NewtonMethods
using Base.Test

# Run tests
tic()
println("Testing NLEQ1")
@time @test include("testNLEQ1.jl")

println("Testing NLEQ2")
@time @test include("testNLEQ2.jl")
toc()
