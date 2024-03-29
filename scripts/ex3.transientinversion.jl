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
newres = 10yr
sum_ind = [a*yr:1yr:a*yr+newres-1yr for a in 1:ustrip(newres):length(τ)-ustrip(newres)+1]
subℳ = cat([sum(ℳ[Ti = At(ind)], dims = Ti) for ind in sum_ind]..., dims = Ti)
arr = cat(Array(reshape(ℳ[At(0yr), :, :], (1, size(ℳ)[2:3]...))), Array(subℳ), dims = 1)
newτ = 0yr:newres:1000yr
ℳnew = DimArray(arr, (Ti(newτ), Modes([m for m in modes]), Cores([c for c in cores])))
ℳ = ℳnew
τ = newτ
res = newres
files = readdir(datadir())
ae_files = [f for f in files if occursin("ae", f)]
d18O_files = [f for f in files if occursin("d18O", f)]
cores = Tuple([Symbol(split(f, ".")[1]) for f in ae_files])

@time e = formatbacon(datadir.(ae_files), datadir.(d18O_files), [c for c in cores], res = res)
#add in measurement uncertainty 
e = fillcovariance!(e, [DiagRule((0.07permil)^2, (:, :))], dims(e.y))

#generate steady-state M

#=
using TMI
TMIversion = "modern_180x90x33_GH11_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
patches = NamedTuple{keys(patches)}([surfacepatch(patches[k]..., γ) for k in keys(patches)])
patches[:East].wet == γ.wet[:, :, 1]
γ.axes
lonlat = (311.0, 61.0)
lonlatindices = (findall(x->x == lonlat[1], γ.lon)[1], findall(x->x == lonlat[2], γ.lat)[1])

tracer = convert(Matrix{Float64}, copy(γ.wet[:, :, 1]))
λ = 5
for i in 1:length(γ.lon)
    for j in 1:length(γ.lat)
        if tracer[i, j] != 0.0 && (i, j) != lonlatindices
            r = sqrt((i - lonlatindices[1])^2 + (j - lonlatindices[2])^2)
            tracer[i, j] = 1-exp(-λ/r)
        end
        
    end
end
bc = BoundaryCondition(tracer, γ.axes[1:2], 0.0, 3, 1, γ.wet[:, :, 1], :name, "longname", "units as string")
patches = NamedTuple{(:East,)}([bc])
M = transientM(corelocs, patches)
=#

σθ = 11K /2 * 2
σδ = 0.8permil / 1 * 2 #std of global surface d18O in WOCE
f = σθ^2/ (20K/permil)
T = [t for t in e.y.dims[1]]
Tᵤ = T[1] - maximum(τ):res:T[end]
u₀ = firstguess(Tᵤ,[m for m in modes], σθ, σδ, f)

#y is a DimArray(T x C)
coeffs = NamedTuple{(:θ, :δ)}([-0.27permil/K, 1])

#for every time, compute
sv = (:θ, :δ)
predict(u) = DimArray(cat([+([+([vec(u[At(t-τi), :, At(s)]' * coeffs[s] * ℳ[At(τi), :, :]) for τi in τ]...) for s in sv]...) for t in T]..., dims = 2)', (Ti(T), Cores([c for c in cores])))

utest = firstguess(Tᵤ, [m for m in modes], σθ, σδ,f, fill_val = [1K, 1permil])
predict(utest.y)

predict(u₀.y)
println("Impulse Response") 
@time E = impulseresponse(predict, u₀.y)

@time yde, ỹde, u₀de, ũ = solvesystem(e, u₀, E, predict);

predict(ũ.x)

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

figure()
scatter(vec(ũ.x[:, :, At(:θ)]), vec(ũ.x[:, :, At(:δ)]), capsize = 5)
lls = linearleastsquares(ustrip.(value.(vec(ũ.x[:, :, At(:δ)]))), ustrip.(value.(vec(ũ.x[:, :, At(:θ)]))))
println("slope [K/permil] = " * string( lls[1]))

figure()
for (i, m) in enumerate(modes)
    subplot(1,5,i)
    for (j, c) in enumerate(cores)
        plot(vec(ℳ[:, At(m), At(c)]), label = c)
    end
end
