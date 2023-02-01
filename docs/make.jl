using Documenter, PhaseGN

makedocs(
    modules  = [PhaseGN],
    sitename = "PhaseGN",
    authors  = "Biplab Mahato",
    pages    = Any[
        "Home"         => "index.md",
        "User Guide"   => "usage.md",
        hide("Indices" => "indices.md")
    ]
)

deploydocs(
    repo = "github.com/biplab37/PhaseGN.git",
)