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
            # "Evaluators" => "forcefield/evaluators.md",
        #     "Page2" => "page2.md",
        #     "Page3" => "page3.md",
        ],
        # "About" => "about.md",
    ]
    #     "Variography" => [
    #         "empirical_variograms.md",
    #         "theoretical_variograms.md",
    #         "fitting_variograms.md"
    #     "Kriging estimators" => "estimators.md",
    #     "Solver comparisons" => "comparisons.md",
    #     "Plotting" => "plotting.md"
    #     ],
    #     "Examples" => "examples.md",
    #     "Contributing" => "contributing.md",
    #     "About" => [
    #     "Community" => "about/community.md",
    #     "License" => "about/license.md",
    #     "Citing" => "about/citing.md"
    #     ]
    # ]
)

