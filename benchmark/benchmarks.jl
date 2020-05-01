# sort out dependencies
using Pkg
tempdir = mktempdir()
Pkg.activate(tempdir)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.add(["BenchmarkTools", "PkgBenchmark"])  # for benchmarking
Pkg.add(["Unitful", "Plots","PlotlyJS", "ORCA"])  # for running experiments
Pkg.resolve()

using BenchmarkTools

const SUITE = BenchmarkGroup()
const PATH_TO_EXAMPLES = "../examples/Epidemiology/"

for file in readdir(joinpath(@__DIR__, PATH_TO_EXAMPLES))
    # temporarily restrict to 2 files only as a prototype
    if file in ["Benchmarking.jl", "Small_SIR.jl"]
        SUITE[file[1:end - length(".jl")]] =
            @benchmarkable include(joinpath(PATH_TO_EXAMPLES, $(file)))
    end
end
