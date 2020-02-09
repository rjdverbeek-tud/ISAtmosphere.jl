using Documenter, ISAtmosphere

makedocs(
    sitename = "ISAtmosphere.jl",
    modules = [ISAtmosphere],
    pages = Any[
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/rjdverbeek-tud/ISAtmosphere.jl",
)
