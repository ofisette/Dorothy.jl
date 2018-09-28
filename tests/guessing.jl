using Test
using Dorothy
using Dorothy: guesselement, guessmass
using Formats
using FormatStreams

datapath = joinpath(@__DIR__, "..", "data")

@testset "Guessing" begin

    m1 = readf("$(datapath)/1BTL.gro")
    guesselements!(m1)
    guessmasses!(m1)
    @test m1.elements[1:4] == ["N", "C", "C", "O"]
    @test m1.masses[1:4] == [14.00, 12.01, 12.01, 15.99]
    @test guessmass("MN1", "") == 0.0
    @test guessmass("HG11", "H") == 1.008

end # @testset
