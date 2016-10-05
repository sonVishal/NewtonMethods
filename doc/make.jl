using Documenter, NewtonMethods

makedocs(modules=NewtonMethods,doctest=true)

deploydocs(deps = Deps.pip("mkdocs","python-markdown-math"),
    repo ="github.com/sonVishal/NewtonMethods.jl.git",
    julia = "0.4.7",
    osname = "linux")
