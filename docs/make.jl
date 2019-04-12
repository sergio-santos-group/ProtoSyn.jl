using Documenter, ProtoSyn
makedocs(
    sitename = "ProtoSyn",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Common" => "common.md",
            "Forcefield" => "forcefield.md",
            "Mutators" => "mutators.md",
            "Drivers" => "drivers.md",
            "Print" => "print.md",
            "Aux" => "aux.md"
        ],
])
deploydocs(repo = "github.com/JosePereiraUA/ProtoSyn.jl.git")
