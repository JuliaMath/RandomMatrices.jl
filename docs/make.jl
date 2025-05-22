push!(LOAD_PATH,"../src/")

using Documenter
using RandomMatrices

DocMeta.setdocmeta!(RandomMatrices, :DocTestSetup, :(using RandomMatrices, Random); recursive=true)

function main()

    makedocs(
    doctest = false,
    clean = true,
    sitename = "RandomMatrices.jl",
    format = Documenter.HTML(
        assets=["assets/init.js"]
    ),
    modules = [RandomMatrices],
    checkdocs = :exports,
    warnonly = false,
    authors = "Andrew Kille",
    pages = [
        "RandomMatrices.jl" => "index.md",
        "Full API" => "API.md"
    ]
    )
end

main()