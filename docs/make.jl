using Documenter, ISAtmosphere

makedocs(
    sitename = "ISAtmosphere.jl",
    modules = [ISAtmosphere],
    pages = Any[
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "rjdverbeek-tud.github.io/ISAtmosphere.jl",
)
