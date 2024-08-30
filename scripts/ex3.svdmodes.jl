#=
Needs to use special environment because ODE needs to be pinned at 6.33.3
which conflicts with some packages I use in the OPTinv.jl environment,
such as DimensionalData 
=#
import Pkg; Pkg.activate("../")
include("../src/OPTinv_alt_env.jl") 

corelocs = core_locations()

surforigin, SVD = generatemodes(corelocs, func = svd)

res = 1 
τ = 0:res:1000
ℳ = transientM(corelocs, SVD.Vt, τ)

jldsave(joinpath("../data/M/svd.jld2"); ℳ, τ, res, surforigin, SVD) #TIME x MODES x CORES
