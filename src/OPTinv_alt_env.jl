#=
Problem: to have stable ODE solutions, I need to add OrdinaryDiffEq to my
environment and pin it to v6.33.3. However, this conflicts with the version
used by DimensionalData, and I really baked DD into OPTinv, so this doesn't work

Workaround for now: run these methods in their own environment to generate the
\scrM matrix 
=#
using CSV, TMItransient, TMI, Statistics, DataFrames,
    Interpolations, Revise, NMF, JLD2, LinearAlgebra

function generatemodes(corelocs; TMIversion = "modern_180x90x33_GH11_GH12", func = svd)
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
    tracers = hcat([surfaceorigin(corelocs[c], Alu, γ).tracer[γ.wet[:, :, 1]] for c in keys(corelocs)]...)
    tracers = 10 .^ tracers
    𝒩 = gaussiandistancematrix(γ, 1, 1000.0)
    mat = tracers' * 𝒩
    if func == nnmf 
        return mat, nnmf(mat, length(corelocs)).H
    elseif func == svd
        return mat,svd(mat)
    else
        return mat,func(mat)
    end    
end

"""
function transientM(corelocs, mat::Matrix, τ; TMIversion = "modern_180x90x33_GH11_GH12")

Propagates user-defined regions (must be TMI-size, 11113) through the L matrix
Returns the globalmean_impulse response observed at core locations 

# Arguments
    - `corelocs`: NamedTuple of core locations as (lat, lon, depth)
    - `mat`: Matrix of surface patterns to propagate through TMItransient L operator
    - `τ`: lag times to compute propagation at

# Output
    - `arr`: 3D Matrix (Ti x cores x modes) 
"""
function transientM(corelocs, mat::Matrix, τ; TMIversion = "modern_180x90x33_GH11_GH12")
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
    
    bcs = [BoundaryCondition(γreshape(mat[i, :],γ), γ.axes[1:2], 0.0, 3, 1, γ.wet[:, :, 1], :i, string(i), string(i)) for i in 1:size(mat)[1]] 
    N = length(corelocs)
    wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,N)
    [wis[i] = interpindex(corelocs[i],γ) for i in 1:N]
    #modes = [k for k in keys(patches)]
    #cores = [c for c in keys(corelocs)]
    arr = Array{Float64}(undef, length(τ), length(bcs), length(corelocs)) # dimensions = Ti x cores x modes
    for (i, bc) in enumerate(bcs)
        #bc = patches[key]
        @time D̄ = stepresponse(TMIversion, bc, γ, L, B, τ, eval_func = observe, args = (wis, γ))
        Ḡ, _ = globalmean_impulseresponse(D̄, τ)
        Ḡ = hcat(Ḡ...)'    
        arr[begin:size(Ḡ)[1], i, :] = Ḡ
    end
    return arr 
end

"""
function transientM(corelocs, patches::NamedTuple, τ; TMIversion = "modern_180x90x33_GH11_GH12")

if you want to use pre-defined TMI regions 
"""
function transientM(corelocs, patches::NamedTuple, τ; TMIversion = "modern_180x90x33_GH11_GH12")
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


function reshape(v::Vector{T}, γ) where T
    #template = zeros(size.(γ.axes)[1][1], size.(γ.axes)[2][1])
    template = Array{T}(undef, size.(γ.axes)[1][1], size.(γ.axes)[2][1])
    template .= NaN
    template[γ.wet[:, :, 1]] .= v
    return template
end


function core_locations()
    locations_file = "../data/EN539_locations.csv"
    cores_df = CSV.read(locations_file,DataFrame)
    locs = [(cores_df.Lon[i], cores_df.Lat[i], cores_df.Depth[i]) for i in 1:length(cores_df.Core)]
    cores = replace.(cores_df.Core, "EN539-"=>"")
    return NamedTuple{Tuple(Symbol.(cores))}(locs)
    
end
