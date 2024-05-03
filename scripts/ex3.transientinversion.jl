import Pkg
Pkg.activate("../")

using OPTinv, Unitful, DimensionalData, BLUEs, Statistics, UnitfulLinearAlgebra, LinearAlgebra, PyPlot, Revise, DrWatson, PyCall, TMI, CSV, DataFrames, JLD2, NaNMath, Dates, RollingFunctions
import Measurements.value as value
import OPTinv.Est


ccrs = pyimport("cartopy.crs")
cm = pyimport("cmocean.cm")
cfeature = pyimport("cartopy.feature")
cu = pyimport("cartopy.util")

#which cores do you want?
corenums_full = [28, 26, 25, 22, 21, 20, 19, 10, 9,13,14]
core_list_full = Symbol.("MC".* string.() .* "A")
#missing  14
corenums = [28, 26, 25, 22, 21, 20, 19, 10, 9,13,14]
#corenums = [28]
#corenums = [26, 22, 13]
#corenums = [26,25]
corenums_sorted = [i for i in corenums_full if i in corenums]
core_list = Symbol.("MC" .* string.(corenums_sorted) .* "A")

#load in ℳ and surface spatial modes 
filename = "nmf.jld2"
ℳnmf, spatialmodesnmf = loadM(filename, core_list)

N = length(corenums_sorted)

#ℳnmf = ℳnmf[:, DimensionalData.Between(1, N), :]
#spatialmodesnmf = spatialmodesnmf[1:N, :]

filename = "svd.jld2"
ℳsvd, spatialmodessvd = loadM(filename, core_list)

y = loadcores(core_list);y.y

# make Cuu matrix
#11K, 0.8permil are global surface std in WOCE 
#σθ = 11K * 2  #* 2.5
#σδ = 1.1permil * 2 #* 2.5     #std of global surface d18O in WOCE
ρ = 0.99
T = Array(y.y.dims[1])
res = unique(diff(T))[1]
Tᵤ = 500.0yr:res:T[end] #arbitrary cut off for Tᵤ, otherwise impulseresponse is very costly
#okay, this is close but I need to work on offdiagrule for a varying σθ
jld = jldopen(datadir("modemags.jld2"))
svdmags, nmfmags = jld["svdmags"], jld["nmfmags"]
σθsvd  = vec(std(svdmags, dims = 1)) * 5K#.* 5K 
σδsvd = σθsvd ./ 10 .* permil/K

σθnmf = vec(std(nmfmags, dims = 1)) * 0.1K  #.* K .* 0.1
σδnmf = σθnmf ./ 10 .* permil/K 


#u₀ =  firstguess(Tᵤ, Array(ℳsvd.dims[2]), σθ, σδ, ρ)
u₀svd = firstguess(Tᵤ, ℳsvd.dims[2][:], σθsvd, σδsvd, ρ)

u₀nmf = firstguess(Tᵤ, ℳnmf.dims[2][:], σθnmf, σδnmf, ρ)
Esvd, predictsvd = loadE("svd.jld2", ℳsvd, Tᵤ, T, σθsvd, σδsvd, ρ,core_list, y.Cnn.ax)
#Enmf, predictnmf = loadE("nmf.jld2", ℳnmf, Tᵤ, T, u₀nmf, core_list, y.Cnn.ax)
Enmf,predictnmf = loadE("nmf.jld2",  ℳnmf, Tᵤ,T, σθnmf, σδnmf, ρ, core_list, y.Cnn.ax)
#Esvd,predictsvd = loadE("svd.jld2",  ℳsvd,Tᵤ,T, σθsvd, σδsvd, ρ, core_list, y.Cnn.ax)

#@time yde, ỹdesvd, u₀de, ũsvd = solvesystem(y, u₀svd, Esvd, predictsvd);
@time yde, ỹdenmf, u₀de, ũnmf = solvesystem(y, u₀nmf, Enmf, predictnmf);
ỹdenmf

fig = figure(figsize = (8, 8))
for (i, core) in enumerate(y.y.dims[2]) 
    subplot(4,3, i)
    plot(yde.x[:, At(core)], label = "y", color = "black")
    #plot(ỹdesvd.x[:, At(core)], label = "ỹ, SVD", color = "red")
    plot(ỹdenmf.x[:, At(core)], label = "ỹ, NMF", color = "blue")
    ylim(-0.25, 0.25)
    title(core)
end
tight_layout()
 #savefig(plotsdir("reconsol.png"))

figure(figsize = (8, 8))
for (i, m) in enumerate(ℳnmf.dims[2])
    for (j, s) in enumerate([:θ, :δ])
        #=
        subplot(2,2,(j-1)*2 + 1)
        if i == 1 
            title("SVD: " * string(s))
        end
        plot(value.(ũsvd.x[:, At(m), At(s)]), label = string(m))
        legend(loc = "center left")
        =#
        subplot(2,2,(j-1)*2 + 2)
        if i == 1 
            title("NMF: " * string(s))
        end        
        plot(value.(ũnmf.x[:, At(m), At(s)]), label = string(m))
        legend(loc = "center left")
    end
end
tight_layout()
savefig(plotsdir("modemags.png"))


#plot just showing correlation between T and δ
#=
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
=#

#plot showing ℳ
#=
figure(figsize = (15,8)) 
for (i, m) in enumerate(ℳ.dims[2])
    subplot(3,4,i)
    for (j, c) in enumerate(y.y.dims[2])
        plot(ℳ[:, At(m), At(c)], label = c)
        xlim([0,200])
        title(m)
    end
end
tight_layout()
=#

TMIversion = "modern_180x90x33_GH11_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
#@time θsvd, δsvd = reconstructsurface(spatialmodessvd, ũsvd)
@time θnmf, δnmf = reconstructsurface(spatialmodesnmf, ũnmf)

#θsvdrs = cat([γreshape(θsvd[i], γ) for i in 1:length(θsvd)]..., dims = 3)
θnmfrs = cat([γreshape(θnmf[i], γ) for i in 1:length(θnmf)]..., dims = 3)
              
#=
#this is reconstructing one time step, takes about 8 seconds, which is very costly for 90 timesteps 
t = 1600.0yr
cd = vec(covariancedims(ũ.dims))
ind = findall(x->x == 1, [x[1] == t && x[3] == :θ for x in cd])
@time sm_big = UnitfulMatrix(Matrix{Float32}(undef, 11113, length(cd)), fill(K, 11113), unit.(vec(ũ.x)));
[sm_big[:, i] = spatialmodes[cd[i][2], :] for i in ind]
@time e = sm_big * ũ


#take a spatial average 
multme = zeros(1, 11113)
multme[1:100] .= 1/100
UnitfulMatrix(multme) * e 

#what if we reconstruct one gridcell and preserve the time covariance?
#this is where the defined data convariance is, so may make more sense

#ULTIMATELY: this is faster but has higer uncertainties for equivalent values 
spatialindex = 1 #which gridcell
multme = zeros(99, 2178)
#multvec = [spatialmodes[i[2], spatialindex] for i in cd]
for (i, t) in enumerate(ũ.dims[1]) #iterate through 99 time indices 
    for (j, c) in enumerate(cd) #iterate through 2178 ũ indices 
        if c[1] == t 
            multme[i, j] = spatialmodes[c[2], spatialindex]
        end 
    end
end
@time test2 = UnitfulMatrix(multme, fill(K, 99), unit.(vec(ũ.x))) * ũ
e.x[1]
test2.x[findall(x->x == 1600.0yr, Array(ũ.dims[1]))[1]]
=#

figure()
clrs = ["red", "blue"]
labels = ["SVD", "NMF"] 
#for (cs, θ, lab) in zip(colors, [θsvdrs, θnmfrs], labels)
for (cs, θ, lab) in zip(clrs, [θnmfrs], labels)
    θ = ustrip.(value.(θ))
    #θ .-= NaNMath.mean(vec(θ))
    θglobal = [NaNMath.mean(vec(θ[:, :, t])) for t in 1:length(Tᵤ)]
    θbox = boxmean(θ, γ.lat, γ.lon, 49, 89, 309, 21)
    plot(DimArray(value.(θglobal .- mean(θglobal)), Ti(Tᵤ)), color = cs, label = lab * ": global")
    plot(DimArray(value.(θbox .- mean(θbox)), Ti(Tᵤ)), color = cs, label = lab* ": 50°W-20°E, 50°N-90°N", linestyle = "dashed")
end


ylabel("Temperature Anomaly from Series Mean [K]")
xlim([ustrip(minimum(T)), 2000])

lmr, lmrt, lmrlon, lmrlat = loadLMR()
lmr .-= NaNMath.mean(lmr)
lmrt = rolling(mean, year.(lmrt), 50)
lmrglobal = rolling(mean, [NaNMath.mean(lmr[:, :, i]) for i in 1:size(lmr)[3]], 50)

lmrbox = rolling(mean, boxmean(lmr, lmrlat, lmrlon, 50, 90, 360-50, 20), 50)
plot(DimArray((lmrglobal .- mean(lmrglobal)) * K, Ti(lmrt)), color = "black", label = "LMR: global", zorder = 0)
plot(DimArray((lmrbox  .- mean(lmrbox))* K, Ti(lmrt)), color = "black", label = "LMR: 50°W-20°E, 50°N-90°N", zorder = 0, linestyle = "dashed" )
legend()
savefig(plotsdir("meants.png"))


res = 5
#plot 50 year means starting at 1500

inds =  findall(x->x == 1000.0yr, Tᵤ)[1]:res:length(Tᵤ)-res
@show length(inds)
lvls = [-0.6:0.1:0.6, -0.6:0.1:0.6]
#stopped plotting delta
savenames = ["theta", "delta"]
#for (var, lev, sn) in zip([θsvd, θnmf], levels, savenames)
for (var, lev, sn) in zip([θnmf], lvls, savenames) 
    figure(figsize = (10, 8))
    anom = var #[v .- var[end] for v in var] 
    for (ii, i) in enumerate(inds)
        ax = subplot(4,5,ii, projection = ccrs.PlateCarree())
        title(string(convert(Int64, ustrip(Tᵤ[i]))) * "-" *string(convert(Int64, ustrip(Tᵤ[i + res]))))
        ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "110m", edgecolor="k", facecolor="gray"))
        ax.set_extent([-80, 30, 0, 90])
        #this is kind of a slow reshape but c'est la vie 
        plotme = γreshape(ustrip.(value.(sum(anom[i:i+res] ./ res))), γ)'
        #add a cyclic point 
        plotme, lon = cu.add_cyclic_point(plotme, coord = ustrip.(γ.lon))
        cf = contourf(lon, γ.lat, plotme, cmap = cm.balance, levels = lev)
        c = contour(lon, γ.lat, plotme, colors = "black", levels = lev)
        ax.clabel(c) 
    end
    tight_layout()
    #savefig(plotsdir(sn * "surfacesol.png"))
end

lab_pt = (360-51, 57)
#Oũsvd = ptobserve(lab_pt, γ, spatialmodessvd, ũsvd) 
Oũnmf = ptobserve(lab_pt, γ, spatialmodesnmf, ũnmf)

df = CSV.read(datadir("TLS1998Fig7Digitized.csv"), DataFrame)
figure(figsize = (16,4))
subplot(1,3,1) 
scatter(df[!, "Salinity"], df[!, "Temperature"], c = df[!, "Year"], cmap = "viridis")
[text(x,y,s) for (x,y,s) in zip(df[!, "Salinity"], df[!, "Temperature"], string.(df[!, "Year"])), vmin = 20, vmax = 95]
plot(df[!, "Salinity"], df[!, "Temperature"], color = "gray") 
xlabel("Salinity")
ylabel("Temperature") 
title("Digitized TLS1998 Figure 7")
#=
subplot(1,3,2)
Stimeseries = vec(Oũsvd.x[DimensionalData.Between(1920yr, 1980yr), At(:δ)]) ./ 0.5255
θtimeseries = vec(Oũsvd.x[DimensionalData.Between(1920yr, 1980yr), At(:θ)])

scatter(value.(Stimeseries), value.(θtimeseries), c = ustrip.(collect(Tᵤ[93:end])), cmap = "viridis", vmin = 1920, vmax = 1995)
[text(x, y, s) for (x,y,s) in zip(ustrip.(value.(Stimeseries)), ustrip.(value.(θtimeseries)), string.(convert.(Int64, ustrip.(Tᵤ[93:end]))))]
xlabel("Salinity anomaly from 1980 [g/kg]")
ylabel("θ anomaly from 1980 [K]")
title("SVD: T, inferred S from 51°W, 57°N")
=#

subplot(1,3,3)
Stimeseries = vec(Oũnmf.x[DimensionalData.Between(1920yr, 1980yr), At(:δ)]) ./ 0.5255
θtimeseries = vec(Oũnmf.x[DimensionalData.Between(1920yr, 1980yr), At(:θ)])
t = collect(1920:10:1980)
scatter(value.(Stimeseries), value.(θtimeseries), c = t, cmap = "viridis", vmin = 1920, vmax = 1995)
[text(x, y, s) for (x,y,s) in zip(ustrip.(value.(Stimeseries)), ustrip.(value.(θtimeseries)), string.(t))]
xlabel("Salinity anomaly from 1980 [g/kg]")
ylabel("θ anomaly from 1980 [K]")
title("NMF: T, inferred S from 51°W, 57°N")
#savefig(plotsdir("labcomp.png"))

figure()
θtimeseries = vec(Oũnmf.x[DimensionalData.Between(1850yr, 1980yr), At(:θ)])
plot(1850:10:1980, ustrip.(value.(θtimeseries)))
meta, data = loadOcean2k()
atl_indices = findall(x->occursin("Atlantic", x), meta[!, "name"])
arc_indices = findall(x->occursin("Arctic", x), meta[!, "name"])
inregion = findall(x->x>50, meta[!, "lat"])

meta[intersect(inregion, vcat(atl_indices, arc_indices)), :]


function convertx(x)
    x[x .== "NaN"] .= NaN
    x = convert(Vector{Float64}, x)
end

        
figure(figsize = (10, 6))
for (ii, i) in enumerate(intersect(inregion, vcat(arc_indices, atl_indices)))
    name = meta[i, "name"] 
    lat = meta[i, "lat"]
    lon = meta[i, "lon"]
    lon = lon < 0 ? lon + 360 : lon
    ptnmf = ptobserve((360 + lon, lat), γ, spatialmodesnmf, ũnmf)
    #ptsvd = ptobserve((lon, lat), γ, spatialmodessvd, ũsvd)
    sst = convertx(data[!, name * "_SST"]) 
    age = convertx(data[!, name * "_Age"])
    subplot(4,2, ii)
    plot(1950 .- age, sst .- NaNMath.mean(sst), color = "black")
    plot(value.(ptnmf.x[:, At(:θ)]), color = "blue")
    #plot(value.(ptsvd.x[:, At(:θ)]), color = "red")
    title(name)
    xlim(ustrip(minimum(T)),2015)
    ylim(-2, 2)
    xlabel("Time [years CE]")
    ylabel("Temperature anomaly [K]")
end
tight_layout()
#savefig(plotsdir("ocean2k.png"))

f = figure()
ax = f.add_subplot(projection = ccrs.PlateCarree())
lats = meta[intersect(inregion, atl_indices), "lat"]
lons = meta[intersect(inregion, atl_indices), "lon"]
names = meta[intersect(inregion, atl_indices), "name"]
ax.coastlines(color = "gray")
ax.set_extent([-30, 30, 50, 70])
scatter(lons, lats)
[text(lon, lat, t[13:end]) for (lon, lat, t) in zip(lons, lats, names)]
#savefig(plotsdir("ocean2klocs.png"))

WHtempnmf = ptobserve((360-45, 57), γ, spatialmodesnmf, ũnmf)
#WHtempsvd = ptobserve((360-45, 57), γ, spatialmodessvd, ũsvd)
EGtempnmf = ptobserve((360-14, 75), γ, spatialmodesnmf, ũnmf)
#EGtempsvd = ptobserve((360-14, 75), γ, spatialmodessvd, ũsvd)

figure();plot(value.(WHtempnmf.x[:, At(:θ)]), color = "red", label = "S. Greenland Anom.")
#plot(value.(WHtempsvd.x[:, At(:θ)]), color = "red", linestyle = "dashed")
xlim(1300,2000)

plot(value.(EGtempnmf.x[:, At(:θ)]), color = "blue", label = "E. Greenland Anom.")
#plot(value.(EGtempsvd.x[:, At(:θ)]), color = "blue", linestyle = "dashed")
legend()
ylim(-1, 1)
twin = gca().twinx()
thornalley = loadThornalley()

for k in keys(thornalley)
    age = thornalley[k][!, "age [CE]"]
    data = thornalley[k][!, "smooth"]
    twin.plot(age, data .- NaNMath.mean(data), color = "black", alpha = 0.5, zorder = 0)
end
ylabel(L"\bar{SS}" * " [mm]")



#Mann's file format is absolutely unhinged 
#=
mann = CSV.read("/home/brynn/Code/OPTinv.jl/data/allproxyfieldrecon", DataFrame, delim = "/t", header = 0)
longlat
=#
