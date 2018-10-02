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
        "Manual" => [
            "Forcefield" => "forcefield.md",
            "Common" => "common.md"
        ],
    ]

)

