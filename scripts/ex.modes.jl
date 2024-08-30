import Pkg;Pkg.activate("../")
using JLD2, Statistics, TMI, OPTinv, PyPlot

jld1 = jldopen("/home/brynn/Code/OPTinv.jl/data/M/svd.jld2")
V = jld1["modes"]
S = jld1["S"]
exp_var = S.^2 / sum(S.^2)
sum(exp_var)

jld2 = jldopen("/home/brynn/Code/OPTinv.jl/data/modemags.jld2")
Vσ = jld2["svdmags"]
std(Vσ, dims = 1)



TMIversion = "modern_180x90x33_GH11_GH12"
γ = Grid(download_ncfile(TMIversion))

θ = γreshape(V' * Vσ[1, :], γ)
cf = contourf(θ');colorbar(cf)

Vσ[1,:]
V[1,:]

maximum(abs.(V[1,:]))
