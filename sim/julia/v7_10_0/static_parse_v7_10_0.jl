const TARGET = joinpath(@__DIR__, "run_fig7_production_v7_10_0.jl")
Meta.parseall(read(TARGET, String))
println("V7.10.0 Julia syntax passed: $(TARGET)")
