import Pkg

filename = "rectangles_1yr.jld2"
directory = "../data/M"
if !isdir(directory)
    mkdir(directory)
end
    Pkg.activate(".")

filepath = joinpath("../data/M", filename)
res = 1

if !isfile(filepath)
    using JLD2 
    include("../src/OPTinv_alt_env.jl")
    surfaceboxes = [[[300, 320], [50,70]], [[320, 340], [50, 70]], [[340, 360],  [50, 70]], [[300, 360], [70, 90]], [[300,360], [30, 50]]]
    modes = (:East,:Central, :West, :North, :South)
    patches = NamedTuple{modes}(surfaceboxes) 
    corelocs = core_locations()
    τ = 0:res:1000
    M = transientM(corelocs, patches, τ)
    M[abs.(M) .> 1] .= 0
    jldsave(filepath; M, τ, surfaceboxes, modes, patches, corelocs)
else
    using JLD2 
    jld = jldopen(filepath)
    M = jld["M"]
    τ = jld["τ"]
    surfaceboxes = jld["surfaceboxes"]
    modes = jld["modes"]
    patches = jld["patches"]
    corelocs = jld["corelocs"]
    close(jld)
end

Pkg.activate("../")

using OPTinv, Unitful, DimensionalData, BLUEs, Statistics, UnitfulLinearAlgebra, LinearAlgebra, PyPlot, Revise, DrWatson
import Measurements.value as value

close("all")

cores = keys(corelocs)
τ = τ * yr
res = unique(diff(τ))[1]
ℳ = formattransientM(M, τ, [m for m in modes], [c for c in cores])

#M has 1 yr resolution, but we're gonna subsample this
#has to be subsampled by summing
res = 10yr 
ℳ, τ = subsampletransientM(ℳ, res)[1]

#read in data at this resolution 
files = readdir(datadir())
ae_files = [f for f in files if occursin("ae", f)]
d18O_files = [f for f in files if occursin("d18O", f)]
listed_cores = Tuple([Symbol(split(f, ".")[1]) for f in ae_files])
#re-order directory cores so they line up with Cores dim 
sort_indices = [findall(x->x== c, cores)[1] for c in listed_cores]
#subset to the cores that show consistency 
subset_cores = Symbol.("MC" .* string.([26, 25, 20, 19, 10 , 9, 13, 14]) .* "A")
subset = [c ∈ subset_cores for c in cores]

#now feed in those files and get e matrix 
@time e = formatbacon(datadir.(ae_files)[sort_indices][subset], datadir.(d18O_files)[sort_indices][subset], [c for c in subset_cores], res = res)
cores = subset_cores

#add in measurement uncertainty 
e = fillcovariance!(e, [DiagRule((0.07permil)^2, (:, :))], dims(e.y))

# make Cuu matrix 
σθ = 11K / 2 
σδ = 0.8permil / 1  #std of global surface d18O in WOCE 
f = σθ^2/ (10K/permil) 
T = [t for t in e.y.dims[1]]
Tᵤ = 1000yr:res:T[end] #arbitrary cut off for Tᵤ, otherwise impulseresponse is very costly
u₀ = firstguess(Tᵤ,[m for m in modes], σθ, σδ, f)

#y is a DimArray(T x C)
coeffs = NamedTuple{(:θ, :δ)}([-0.27permil/K, 1])

#for every time, compute
sv = (:θ, :δ)
predict(u) = DimArray(cat([+([+([vec(u[At(t-τi), :, At(s)]' * coeffs[s] * ℳ[At(τi), :, At(cores)]) for τi in τ[t .- τ .> minimum(Tᵤ)]]...) for s in sv]...) for t in T]..., dims = 2)', (Ti(T), Cores([c for c in cores])))

println("Impulse Response")
filename = split(filename, ".")[1] .* "_E.jld2"
filepath = joinpath("../data/M", filename) 
if !isfile(filepath) 
    @time E = impulseresponse(predict, u₀.y) #328 seconds
    jldsave(filepath; E) 
else
    println("Loading in pre-saved file") 
    jld = jldopen(filepath)
    E = jld["E"]
    close(jld) 
end 

@time yde, ỹde, u₀de, ũ = solvesystem(e, u₀, E, predict);

fig = figure(figsize = (8, 8))
for (i, core) in enumerate(cores) 
    subplot(4,3, i)
    plot(yde.x[:, At(core)], label = "y", color = "black")
    plot(ỹde.x[:, At(core)], label = "ỹ", color = "red")
    ylim(-0.25, 0.25)
    title(core)
end
tight_layout()
 
figure(figsize = (8,4))
for (i, m) in enumerate(modes)
    for (j, s) in enumerate([:θ, :δ])
        subplot(2, length(modes), (j - 1) * length(modes) + i)
        title(string(m) * " " * string(s)) 
        plot(u₀de.x[:, At(m), At(s)], label = "u₀", color = "gray")
        plot(ũ.x[:, At(m), At(s)], label = "ũ", color = "red")    
    end
end
tight_layout()

figure(figsize = (10, 5))
subplot(1,2,1) 
scatter(vec(ũ.x[:, :, At(:θ)]), vec(ũ.x[:, :, At(:δ)]), capsize = 5)
subplot(1,2,2)
halfway = convert(Int64, length(ũ.v)/2)
scatter(ustrip.(vec(ũ.v)[begin:halfway]), ustrip.(vec(ũ.v)[halfway+1:end]))
lls = linearleastsquares(ustrip.(value.(vec(ũ.x[:, :, At(:δ)]))), ustrip.(value.(vec(ũ.x[:, :, At(:θ)]))))
println("slope [K/permil] = " * string(lls[1]))

figure(figsize = (10,3)) 
for (i, m) in enumerate(modes)
    subplot(1,5,i)
    for (j, c) in enumerate(cores)
        plot(ustrip.(τ), vec(ℳ[:, At(m), At(c)]), label = c)
        title(m)
    end
end
tight_layout()



