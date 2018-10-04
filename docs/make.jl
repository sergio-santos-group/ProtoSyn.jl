if Base.HOME_PROJECT[] !== nothing
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
  end

using Documenter, ProtoSyn

push!(LOAD_PATH,"../src/")

makedocs(
    format = :html,
    sitename = "ProtoSyn.jl",
    pages = [
        "Home" => "index.md",
        "Guide" => "guide.md",
        "Manual" => [
            "Common" => "common.md",
            "Forcefield" => "forcefield.md",
            "Mutators" => "mutators.md",
            "Drivers" => "drivers.md",
            "Print" => "print.md",
            "Aux" => "aux.md",
            "Input JSON Schemas" => "json.md"
        ],
    ]
)

deploydocs(
    repo = "github.com/sergio-santos-group/ProtoSyn.jl.git",
    osname = "linux",
    julia = "1.0",
    deps = nothing,
    make = nothing,
    target = "build",
)
