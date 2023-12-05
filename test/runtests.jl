using VaidyaPT
using Test

Ψ = Field(parfile = "test/Schwarzschild_pars.txt")     
Initialize(Ψ)                                             
Evolve(Ψ, "test/test_data/") 
data = read_data("test/test_data/")
    
@testset "VaidyaPT.jl" begin
    @test data.v[end] == 50.0
    @test abs(data.ψ[end] - 0.03195) < 0.005              
end
