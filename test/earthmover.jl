using OptimalTransport
using Distances

import LinearAlgebra
import Random
using Test

Random.seed!(1234)

@testset "Initial map" begin
    # example from "Computational Optimal Transport"
    P = zeros(3, 3)
    OptimalTransport.initialmap!(P, [0.2, 0.5, 0.3], [0.5, 0.1, 0.4])
    @test P == [0.2 0 0; 0.3 0.1 0.1; 0 0 0.3]
end

# example taken from https://www.me.utexas.edu/~jensen/models/network/net8.html
@testset "Simple problem" begin
    C = [3 1 1000; 4 2 4; 1000 3 3]
    P, f, g = earthmover([5, 7, 3], [7, 3, 5], C)

    @test vec(sum(P; dims=1)) == [7, 3, 5]
    @test vec(sum(P; dims=2)) == [5, 7, 3]
    @test LinearAlgebra.dot(C, P) == 46
end

@testset "Total variation" begin
    histogram1 = LinearAlgebra.normalize!(rand(20), 1)
    histogram2 = LinearAlgebra.normalize!(rand(20), 1)

    C = (!==).(1:20, (1:20)')
    P, f, g = earthmover(histogram1, histogram2, C)

    @test vec(sum(P; dims=1)) ≈ histogram2
    @test vec(sum(P; dims=2)) ≈ histogram1
    @test LinearAlgebra.dot(C, P) ≈ totalvariation(histogram1, histogram2)
end