push!(LOAD_PATH, "../")

using Documenter
using Simulation

makedocs(
    modules = [Simulation],
    sitename = "Simulation.jl",
    format=Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages = [
        "Home" => "index.md",
        "Design" => "Design.md",
        "Structure" => "Structure.md",
        "Structure_gh" => "Structure_gh.md",
    ],
    strict=true,
    checkdocs=:none,
)
