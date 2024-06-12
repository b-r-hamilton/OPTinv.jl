using OPTinv
using Test

@testset "OPTinv.jl" begin

    @testset "Data Load" begin
        lab_pt = (360-51, 57)
        reg_box = [[309, 21], [49, 89]]
        en4time, en4temp, en4salinity = loadEN4("/home/brynn/Code/oceanFP/data/EN4/analyses",  lab_pt...)
        en4time, en4temp, en4salinity = loadEN4("/home/brynn/Code/oceanFP/data/EN4/analyses",  reg_box...)

    end
    
end
