using Simulation
using Test
using Distributions
using Unitful.DefaultSymbols
using Simulation.Units

include("TestCases.jl")

epi = TestEpiSystem()

@test sum(epi.abundances.matrix, dims = 2)[:, 1] == epi.epilist.abun
@test_nowarn gettraitrel(epi)
@test gettraitrel(epi) == epi.relationship
@test_nowarn gethabitat(epi)
@test gethabitat(epi) == epi.epienv.habitat
@test_nowarn getsize(epi)
@test getsize(epi) == size(epi.abundances.matrix, 2) .* epi.epienv.habitat.size^2
@test_nowarn getgridsize(epi)
@test getgridsize(epi) == epi.epienv.habitat.size
@test_nowarn getdispersaldist(epi, 1)
@test getdispersaldist(epi, 1) == epi.epilist.movement.kernels[1].dist
@test_nowarn getdispersaldist(epi, "Virus")
@test getdispersaldist(epi, "Virus") == epi.epilist.movement.kernels[1].dist
@test_nowarn getdispersalvar(epi, 1)
@test_nowarn getdispersalvar(epi, "Virus")

@testset "Comparison" begin
    @testset "Equality" begin
        epi_2 = deepcopy(epi)
        @test epi_2 == epi
    end
    @testset "Approximate equality" begin
        epi_2 = deepcopy(epi)
        @test isapprox(epi_2, epi)
        # Make some small change within the approximate equality range
        epi_2.abundances.matrix[2, 1] += 1
        @test epi_2 != epi
        @test isapprox(epi_2, epi)
        # Make some larger change outside the approximate equality range
        epi_2.abundances.matrix[2, 1] += 1000
        @test epi_2 != epi
        @test !isapprox(epi_2, epi)
    end
end
