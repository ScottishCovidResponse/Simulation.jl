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
        "Model Structure" => "model_structure.md",
        "Model Development" => "model_development.md",
        "API" => "api.md",
    ],
    strict=true,
    checkdocs=:none,
)
