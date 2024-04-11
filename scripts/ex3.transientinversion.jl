
import Pkg
Pkg.activate("../")

using OPTinv, Unitful, DimensionalData, BLUEs, Statistics, UnitfulLinearAlgebra, LinearAlgebra, PyPlot, Revise, DrWatson, PyCall, JLD2
import Measurements.value as value
using DimensionalData: @dim
@dim Gridcells "Grid cells"
ccrs = pyimport("cartopy.crs")
cm = pyimport("cmocean.cm")
cfeature = pyimport("cartopy.feature")
cu = pyimport("cartopy.util")

filename = "arr.jld2"
directory = "../data/M"
filepath = joinpath("../data/M", filename)

jld = jldopen(filepath)
M = jld["arr"]
τ = jld["τ"]
modes = 1:11
cores = keys(core_locations())
close(jld)

τ = τ * yr
res = unique(diff(τ))[1]
ℳ = formattransientM(M, τ, [m for m in modes], [c for c in cores])

#M has 1 yr resolution, but we're gonna subsample this
#has to be subsampled by summing
res = 10yr 
ℳ, τ = subsampletransientM(ℳ, res)

#read in data at this resolution 
files = readdir(datadir())
ae_files = [f for f in files if occursin("ae", f)]
d18O_files = [f for f in files if occursin("d18O", f)]
listed_cores = Tuple([Symbol(split(f, ".")[1]) for f in ae_files])
#re-order directory cores so they line up with Cores dim 
sort_indices = [findall(x->x== c, cores)[1] for c in listed_cores]
#=
#subset to the cores that show consistency 
subset_cores = Symbol.("MC" .* string.([26, 25, 20, 19, 10 , 9, 13, 14]) .* "A")
subset_cores = Symbol.("MC" .* string.([28, 10, 14]) .* "A")
subset = [c ∈ subset_cores for c in cores]

#now feed in those files and get e matrix 
@time e = formatbacon(datadir.(ae_files)[sort_indices][subset], datadir.(d18O_files)[sort_indices][subset], [c for c in subset_cores], res = res)
cores = subset_cores
=#
@time e = formatbacon(datadir.(ae_files)[sort_indices], datadir.(d18O_files)[sort_indices], [c for c in cores], res = res)
#add in measurement uncertainty 
e = fillcovariance!(e, [DiagRule((0.07permil)^2, (:, :))], dims(e.y))

# make Cuu matrix
#11K, 0.8permil are global surface std in WOCE 
σθ = 11K /2 
σδ = 1.1permil /2     #std of global surface d18O in WOCE
ρ = 0.99
T = [t for t in e.y.dims[1]]
Tᵤ = 1000yr:res:T[end] #arbitrary cut off for Tᵤ, otherwise impulseresponse is very costly
u₀ = firstguess(Tᵤ,[m for m in modes], σθ, σδ, ρ)

#y is a DimArray(T x C)
coeffs = NamedTuple{(:θ, :δ)}([-0.27permil/K, 1])

#for every time, compute
sv = (:θ, :δ)
#might need At(cores) for subsetting ℳ
predict(u) = DimArray(cat([+([+([vec(u[At(t-τi), :, At(s)]' * coeffs[s] * ℳ[At(τi), :, :]) for τi in τ[t .- τ .> minimum(Tᵤ)]]...) for s in sv]...) for t in T]..., dims = 2)', (Ti(T), Cores([c for c in cores])))
vec(predict(u₀.y))
vec(e.y)

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
 
figure(figsize = (16,4))
for (i, m) in enumerate(modes)
    for (j, s) in enumerate([:θ, :δ])
        subplot(2, length(modes), (j - 1) * length(modes) + i)
        title(string(m) * " " * string(s)) 
        plot(u₀de.x[:, At(m), At(s)], label = "u₀", color = "gray")
        plot(ũ.x[:, At(m), At(s)], label = "ũ", color = "red")    
    end
end
tight_layout()

figure(figsize = (20, 5))
subplot(1,2,1) 
scatter(vec(ũ.x[:, :, At(:θ)]), vec(ũ.x[:, :, At(:δ)]), capsize = 5)
subplot(1,2,2)
halfway = convert(Int64, length(ũ.v)/2)
scatter(ustrip.(vec(ũ.v)[halfway+1:end]), ustrip.(vec(ũ.v)[begin:halfway]))
ylabel("Variations in Temperature [K]")
xlabel("Variations in d18O [‰]") 
lls = linearleastsquares(ustrip.(value.(vec(ũ.x[:, :, At(:δ)]))), ustrip.(value.(vec(ũ.x[:, :, At(:θ)]))))
println("slope [K/permil] = " * string(lls[1]))

figure(figsize = (8,8)) 
for (i, m) in enumerate(modes)
    subplot(3,4,i)
    for (j, c) in enumerate(cores)
        plot(ustrip.(τ), vec(ℳ[:, At(m), At(c)]), label = c)
        title(m)
    end
end
tight_layout()

Pkg.activate(".")
include("../src/OPTinv_alt_env.jl")
mat = generatemodes(core_locations())
Pkg.activate("../")

TMIversion = "modern_180x90x33_GH11_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

mat' * vec(ũ.x[At(1600.0yr), :, At(:θ)])

θvals = [mat' * vec(ũ.x[At(t), :, At(:θ)]) for t in Tᵤ]
δvals = [mat' * vec(ũ.x[At(t), :, At(:δ)]) for t in Tᵤ]
#θ = cat([reshape(mat' * vec(ũ.x[At(t), :, At(:θ)]), γ) for t in Tᵤ]..., dims = 3)
#δ = cat([reshape(mat' * vec(ũ.x[At(t), :, At(:δ)]), γ) for t in Tᵤ]..., dims = 3)

θmean, θmin, θmax = [[f(θvals[t]) for t in 1:length(Tᵤ)] for f in [mean, minimum, maximum]]
δmean, δmin, δmax  = [[f(δvals[t]) for t in 1:length(Tᵤ)] for f in [mean, minimum, maximum]]

figure()
subplot(2,1,1)
#fill_between(ustrip.(Tᵤ), y1 = θmean .- θstd, y2 = θmean .+ θstd, color = "blue", alpha = 0.5)
plot(DimArray(θmean, Ti(Tᵤ)), color = "blue")
ylabel("T [K]")
xlim([1400, 2000])
subplot(2,1,2)
plot(DimArray(δmean, Ti(Tᵤ)), color = "blue")
ylabel("δ¹⁸O [‰]")
xlabel("Time [years CE]")
xlim([1400, 2000])

lev = 5
res = 5
inds =  51:res:length(Tᵤ)-res
levels = [-0.6:0.1:0.6, -0.06:0.01:0.06] 
for (var, lev) in zip([θvals, δvals], levels) 
    figure(figsize = (10, 8))
    anom = [v .- var[end] for v in var] 
    for (ii, i) in enumerate(inds)
        ax = subplot(3,3,ii, projection = ccrs.PlateCarree())
        title(string(convert(Int64, ustrip(Tᵤ[i]))) * "-" *string(convert(Int64, ustrip(Tᵤ[i + res]))))
        #ax.coastlines()
        ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "110m", edgecolor="k", facecolor="gray"))
        ax.set_extent([-80, 30, 0, 90])
        plotmevals = mean(hcat(anom[i:i+res]...), dims = 2)[:,1] #vector
        plotme = reshape(ustrip.(value.(plotmevals)), γ)'
        display(plotme)
        #plotme, lon = cu.add_cyclic_point(plotme, coord = ustrip.(γ.lon))
        lon = γ.lon
        cf = contourf(lon, γ.lat, plotme, cmap = cm.balance, levels = lev)#, vmin = -vm, vmax = vm, levels = lev)
        c = contour(lon, γ.lat, plotme, colors = "black", levels = lev)#, vmin = -vm, vmax = vm, levels = lev )
        ax.clabel(c) 
        if i == inds[end]
            #colorbar(cf)
        end
    end
    tight_layout()
end


surface = convert.(Int64, zeros(180,90))
surface[γ.wet[:, :, 1]] .= 1:11113

lab_pt = [360-51, 57]
lab_spot = [findall(x->x == lab_pt[1],γ.lon)[1],findall(x->x==lab_pt[2],γ.lat)[1]]
ind = surface[lab_spot[1], lab_spot[2]]
mapda = DimArray(mat[:, ind], (Modes(1:11)))



stmp = DimArray(zeros(99, 2), (Ti(Tᵤ), StateVar([:θ, :δ])))
x = vec(covariancedims(stmp.dims))
y = vec(covariancedims(ũ.x.dims))
tmat = zeros(length(stmp), length(ũ.x))
for (i, xi) in enumerate(x)
    for (j, yj) in enumerate(y)
        if xi[1] == yj[1]
            if xi[2] == yj[3]
                tmat[i,j] = mapda[yj[2]]
            end
        end
    end
end

obs = UnitfulMatrix(tmat, vcat(fill(K, 99), fill(permil, 99)), unitrange(ũ.v)) * ũ
obsde = DimEstimate(obs.v, obs.C, (Ti(Tᵤ), StateVar([:θ, :δ])))

df = CSV.read(datadir("TLS1998Fig7Digitized.csv"), DataFrame)
figure(figsize = (9,4))
subplot(1,2,1) 
scatter(df[!, "Salinity"], df[!, "Temperature"], c = df[!, "Year"], cmap = "viridis")
[text(x,y,s) for (x,y,s) in zip(df[!, "Salinity"], df[!, "Temperature"], string.(df[!, "Year"])), vmin = 20, vmax = 95]
plot(df[!, "Salinity"], df[!, "Temperature"], color = "gray") 
xlabel("Salinity")
ylabel("Temperature") 
title("Digitized TLS1998 Figure 7") 
subplot(1,2,2)
Stimeseries = vec(obsde.x[DimensionalData.Between(1920yr, 1980yr), At(:δ)]) ./ 0.5255
θtimeseries = vec(obsde.x[DimensionalData.Between(1920yr, 1980yr), At(:θ)])
θtimeseries =[θvals[i][ind] for i in 93:length(θvals)]
#s = scatter(Stimeseries, θtimeseries)
scatter(value.(ustrip.(Stimeseries)), value.(ustrip.(θtimeseries)), c = ustrip.(collect(Tᵤ[93:end])), cmap = "viridis", vmin = 1920, vmax = 1980)
[text(x, y, s) for (x,y,s) in zip(ustrip.(value.(Stimeseries)), ustrip.(value.(θtimeseries)), string.(convert.(Int64, ustrip.(Tᵤ[93:end]))))]
xlabel("Salinity anomaly from 1980 [g/kg]")
ylabel("θ anomaly from 1980 [K]")
title("T, inferred S from 51°W, 57°N")
