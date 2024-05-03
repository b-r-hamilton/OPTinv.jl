#=
Needs to use special environment because ODE needs to be pinned at 6.33.3
which conflicts with some packages I use in the OPTinv.jl environment,
such as DimensionalData 
=#
import Pkg; Pkg.activate(".")
include("../src/OPTinv_alt_env.jl") 

corelocs = core_locations()

modes = generatemodes(corelocs, func = nnmf)
res = 1 
τ = 0:res:1000
arr = transientM(corelocs, modes, τ)

jldsave(joinpath("../data/M/svd.jld2"); arr, τ, modes) #TIME x MODES x CORES

figure()
for i in 1:11
    subplot(3,4,i)
    for j in 1:11
        plot(arr[:, j, i])
    end
    xlim([-5,200])
end
