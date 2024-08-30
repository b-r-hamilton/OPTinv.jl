import Pkg
Pkg.activate("../")

using OPTinv, Unitful, DimensionalData, BLUEs, Statistics, UnitfulLinearAlgebra, LinearAlgebra, PyPlot, Revise, DrWatson, PyCall, TMI, CSV, DataFrames, JLD2, NaNMath, Dates, RollingFunctions, PaleoData, DateFormats, Measurements, GH19
import Measurements.value as value
import OPTinv.Est

#plotting python packages 
ccrs = pyimport("cartopy.crs")
cm = pyimport("cmocean.cm")
cfeature = pyimport("cartopy.feature")
cu = pyimport("cartopy.util")
mpath = pyimport("matplotlib.path")
patches = pyimport("matplotlib.patches")
cmap = pyimport("matplotlib.cm")

allcores = Symbol.("MC" .* string.([28, 26, 25, 22, 21, 20, 19, 10, 9,13,14]) .* "A")
oldcores = Symbol.("MC" .* string.([26, 22, 13]) .* "A")

#= We want to put a statistical constraint that temperature should decrease regionally
according to my Oc2k estimate (0.36 ± 0.13 °C) 
=#


#these scraps are scattered throughout OPTinv.invert
x = inversion(oldcores, "svd", "old", "blue")
corenums_full = Symbol.("MC" .* string.([28, 26, 25, 22, 21, 20, 19, 10, 9,13,14]) .* "A")
corenums_sorted = [i for i in corenums_full if i in x.cores]
filename = x.modetype * ".jld2"
ℳ, spatialmodes = loadM(filename, corenums_sorted, normalize = false)
y = loadcores(corenums_sorted)
ρ = 0.99
T = Array(y.y.dims[1])
res = unique(diff(T))[1]

TMIversion = "modern_180x90x33_GH11_GH12"
γ = Grid(download_ncfile(TMIversion))
Tᵤ = 500.0yr:res:T[end]



#=
if x.modetype == "svd"    
σθ *= 5
else
σθ *= 0.1
end
=#

#σθ *= 5
#σθ = vec(2K ./ jldopen(datadir("M/" * x.modetype * ".jld2"))["absmax"])

# =========== MODIFIED `firstguess` CODE ======= #
modes = 1:11
L = sum(γ.wet[:, :, 1])
V = UnitfulMatrix(spatialmodes, fill(K, size(spatialmodes)[1]), fill(K, size(spatialmodes)[2])) #VVᵀ = I
Tmax = -500yr:10yr:2020yr
oc2kcool = collect((-0.37±0.26)/1000 * K/yr .* Tmax)
sol_anom_ind = findall(x->1850yr<x<1980yr, Tmax)
oc2kcool .-= mean(oc2kcool[sol_anom_ind])
plot(DimArray(oc2kcool, Ti(Tmax)))


jld = jldopen(OPTinv.datadir("modemags.jld2")) #NOT NORMALIZED
mags = jld[x.modetype * "mags"]
σθ = vec(std(mags, dims = 1)) * K

σθm = UnitfulMatrix(ustrip.(diagm(σθ.^2)), fill(K, 11), fill(K^-1, 11))
σδ =  σθ ./ 10 .* permil/K
nmat = Matrix{Any}(undef, length(Tmax), length(modes))
for (i, t) in enumerate(Tmax)
    m = UnitfulMatrix(fill(value(oc2kcool[i]), L))
    Cmm = UnitfulMatrix(Diagonal(fill(ustrip(Measurements.uncertainty(oc2kcool[i])^2), L)), fill(K, L), fill(K^-1, L))
    n = m'*transpose(V)
    Cnn = V*Cmm*transpose(V) + σθm
    if i == 1 display(Cnn) end
    nmat[i, :] = measurement.(vec(n), sqrt.(diag(parent(Cnn)))K)
end
da = DimArray(nmat, (Ti(Tmax), Modes(modes)))

u₀ = firstguess(Tᵤ, ℳ.dims[2][:], da, σδ, ρ) #u₀ with correct T_u
E, predict = loadE(filename, ℳ, Tᵤ, T, da, σδ, ρ, corenums_sorted, y.Cnn.ax)
yde, ỹde, u₀de, ũ = solvesystem(y, u₀, E, predict)
#θ, δ = reconstructsurface(spatialmodes, ũ)

ỹ₀ = predict(u₀de.x)

figure(figsize = (8,4))
for (i,c) in enumerate(oldcores)
    subplot(1,3,i)
    plot(yde.x[:, At(c)], color = "black")
    plot(ỹde.x[:, At(c)], color = "blue")
    title(c)
end

# ============================= REGION MEAN ===================== # 
regionmeanindices = γbox(γ, 49, 89, 309, 21)
figure(figsize = (8,4))
subplot(121)
θmean = estimate(ũ, spatialmodes, :θ, spatialinds = regionmeanindices) 
plot(θmean)
ylim(-2.5, 2.5)
title("Regional Mean Temperature")
subplot(122)
δmean = estimate(ũ, spatialmodes, :δ, spatialinds = regionmeanindices)
plot(δmean)
title("Regional Mean d18O")
ylim(-0.3, 0.3)
figure()
s = scatter(Array(δmean), Array(θmean), c = ustrip.(collect(Tᵤ)))
xl = gca().get_xlim()
x = range(xl[1] + 0.03, xl[2]-0.03, length = 100)
yl = gca().get_ylim()
plot(xl, xl .* 10)

#colorbar(s) #broken for some reason? 

# =================== MAPS ================== #
