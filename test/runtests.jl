using OPTinv
using Test

@testset "OPTinv.jl" begin
    @testset "offsets" begin
        include("../scripts/offsets.jl")
    end
end
