using Pkg


# Can't include this in Project.toml because for some Reaso
# get the error expected SimulationData to be registered
# https://discourse.julialang.org/t/cant-instantiate-project-with-two-unregistered-packages/36703/5


println("Adding SimulationData...")
Pkg.add(PackageSpec(url="https://github.com/ScottishCovidResponse/SimulationData.jl"))
Pkg.build("SimulationData")
