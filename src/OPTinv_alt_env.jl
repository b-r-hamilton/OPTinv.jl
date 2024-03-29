#=
Problem: to have stable ODE solutions, I need to add OrdinaryDiffEq to my
environment and pin it to v6.33.3. However, this conflicts with the version
used by DimensionalData, and I really baked DD into OPTinv, so this doesn't work

Workaround for now: run these methods in their own environment to generate the
\scrM matrix 
=#
using CSV, TMItransient, TMI, Statistics, DataFrames, Interpolations, Revise

"""
copy and pasted from OPTinv.jl (yes I know its bad practice) 
"""
function core_locations()
    
    locations_file = "../data/EN539_locations.csv"
    cores_df = CSV.read(locations_file,DataFrame)
    locs = [(cores_df.Lon[i], cores_df.Lat[i], cores_df.Depth[i]) for i in 1:length(cores_df.Core)]
    cores = replace.(cores_df.Core, "EN539-"=>"")
    return NamedTuple{Tuple(Symbol.(cores))}(locs)
    
end

function transientM(corelocs, patches, τ; TMIversion = "modern_180x90x33_GH11_GH12")
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
    patches = NamedTuple{keys(patches)}([surfacepatch(patches[k]..., γ) for k in keys(patches)])
    N = length(corelocs)
    wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,N)
    [wis[i] = interpindex(corelocs[i],γ) for i in 1:N]
    modes = [k for k in keys(patches)]
    cores = [c for c in keys(corelocs)]
    arr = Array{Float64}(undef, length(τ), length(modes), length(cores)) # dimensions = Ti x cores x modes
    for (i, key) in enumerate(keys(patches))
        bc = patches[key]
        @time D̄ = stepresponse(TMIversion, bc, γ, L, B, τ, eval_func = observe, args = (wis, γ))
        Ḡ, _ = globalmean_impulseresponse(D̄, τ)
        Ḡ = hcat(Ḡ...)'    
        arr[begin:size(Ḡ)[1], i, :] = Ḡ
    end
    #ℳ = DimArray(arr, (Ti(τ * yr), Modes(modes), Cores(cores)))
    #return ℳ
    return arr 
end
