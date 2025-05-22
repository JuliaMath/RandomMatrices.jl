using Documenter
using RandomMatrices
using Test

@testset "Doctests" begin
    DocMeta.setdocmeta!(RandomMatrices, :DocTestSetup, :(using RandomMatrices, Random); recursive=true)
    doctest(RandomMatrices)
end