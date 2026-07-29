using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
const REQUIRED = [
    "sim/julia/Project.toml",
    "sim/julia/Manifest.toml",
    "sim/julia/fig7_geometry_contract.toml",
    "sim/julia/v5/run_fig7_fluid_domain_area_mismatch_v5.jl",
    "sim/julia/v7_7_0/run_fig7_local_leaflet_interface_v7_7_0.jl",
    "sim/julia/v7_10_0/run_fig7_production_v7_10_0.jl",
]

for relative_path in REQUIRED
    path = joinpath(ROOT, relative_path)
    isfile(path) || error("Required file is missing: $relative_path")
end

for relative_path in filter(path -> endswith(path, ".jl"), REQUIRED)
    Meta.parseall(read(joinpath(ROOT, relative_path), String))
end

contract = TOML.parsefile(joinpath(
    ROOT, "sim", "julia", "fig7_geometry_contract.toml"))
isempty(contract) && error("The Fig. 7 geometry contract is empty.")

println("Static Julia dependency and syntax checks passed.")
