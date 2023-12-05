using VaidyaPT
using Test

cd(@__DIR__)

Ψ = Field(parfile = "Schwarzschild_pars.txt")     
Initialize(Ψ)                                             
Evolve(Ψ, "test_data/") 
data = read_data("test_data/")
    
@testset "VaidyaPT.jl" begin
    @test data.v[end] == 50.0
    @test abs(data.ψ[end] - 0.03195) < 0.005              
end
