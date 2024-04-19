
import Pkg
Pkg.activate("../")

using OPTinv, Unitful, DimensionalData, BLUEs, Statistics, UnitfulLinearAlgebra, LinearAlgebra, PyPlot, Revise, DrWatson, PyCall, TMI, CSV, DataFrames, JLD2
import Measurements.value as value
import OPTinv.Est


ccrs = pyimport("cartopy.crs")
cm = pyimport("cmocean.cm")
cfeature = pyimport("cartopy.feature")
cu = pyimport("cartopy.util")

#which cores do you want? 
core_list = Symbol.("MC".* string.([28, 26, 25, 22, 21, 20, 19, 10, 9,13, 14]) .* "A")
#load in ℳ and surface spatial modes 
filename = "nmf.jld2"
ℳnmf, spatialmodesnmf = loadM(filename, core_list)
filename = "svd.jld2"
ℳsvd, spatialmodessvd = loadM(filename, core_list)

y = loadcores(core_list)

# make Cuu matrix
#11K, 0.8permil are global surface std in WOCE 
σθ = 11K * 2  #* 2.5
σδ = 1.1permil * 2 #* 2.5     #std of global surface d18O in WOCE
ρ = 0.99
T = Array(y.y.dims[1])
res = unique(diff(T))[1]
Tᵤ = 1000yr:res:T[end] #arbitrary cut off for Tᵤ, otherwise impulseresponse is very costly

u₀svd = firstguess(Tᵤ, Array(ℳsvd.dims[2]), σθ, σδ, ρ)
u₀nmf = firstguess(Tᵤ, Array(ℳsvd.dims[2]), σθ*20, σδ*20, ρ)
Esvd, predictsvd = loadE("svd.jld2", ℳsvd, Tᵤ, T, u₀svd, core_list, y.Cnn.ax)
Enmf, predictnmf = loadE("nmf.jld2", ℳnmf, Tᵤ, T, u₀nmf, core_list, y.Cnn.ax) 
@time yde, ỹdesvd, u₀de, ũsvd = solvesystem(y, u₀svd, Esvd, predictsvd);
@time _, ỹdenmf, _, ũnmf = solvesystem(y, u₀nmf, Enmf, predictnmf);


fig = figure(figsize = (8, 8))
for (i, core) in enumerate(y.y.dims[2]) 
    subplot(4,3, i)
    plot(yde.x[:, At(core)], label = "y", color = "black")
    plot(ỹdesvd.x[:, At(core)], label = "ỹ, SVD", color = "red")
    plot(ỹdenmf.x[:, At(core)], label = "ỹ, NMF", color = "blue")
    ylim(-0.25, 0.25)
    title(core)
end
tight_layout()

figure(figsize = (8, 8))
for (i, m) in enumerate(ℳsvd.dims[2])
    for (j, s) in enumerate([:θ, :δ])
        subplot(2,2,(j-1)*2 + 1)
        if i == 1 
            title("SVD: " * string(s))
        end
        plot(value.(ũsvd.x[:, At(m), At(s)]), label = string(m))
        legend(loc = "center left")
        subplot(2,2,(j-1)*2 + 2)
        if i == 1 
            title("NMF: " * string(s))
        end        
        plot(value.(ũnmf.x[:, At(m), At(s)]), label = string(m))
        legend(loc = "center left")
    end
end
tight_layout()

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
@time θsvd, δsvd = reconstructsurface(spatialmodessvd, ũsvd)
@time θnmf, δnmf = reconstructsurface(spatialmodesnmf, ũnmf)

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

#=
θμbox, δμbox = boxmean(θ, δ, γ, 49, 89, 309, 21)

θmean, θmin, θmax = [[f(θ[t]) for t in 1:length(Tᵤ)] for f in [mean, minimum, maximum]]
δmean, δmin, δmax  = [[f(δ[t]) for t in 1:length(Tᵤ)] for f in [mean, minimum, maximum]]

figure()
subplot(2,1,1)
plot(DimArray(value.(θmean), Ti(Tᵤ)), color = "blue", label = "global")
plot(DimArray(value.(θμbox), Ti(Tᵤ)), color = "red", label = "50°W-20°E, 50°N-90°N")
legend()
ylabel("T [K]")
xlim([1400, 2000])
subplot(2,1,2)
plot(DimArray(value.(δmean), Ti(Tᵤ)), color = "blue", label = "global")
plot(DimArray(value.(δμbox), Ti(Tᵤ)), color = "red", label = "50°W-20°E, 50°N-90°N")
ylabel("δ¹⁸O [‰]")
xlabel("Time [years CE]")
xlim([1400, 2000])
=#

res = 5
#plot 50 year means starting at 1500 
inds =  51:res:length(Tᵤ)-res
levels = [-0.6:0.1:0.6, -0.6:0.1:0.6]
#stopped plotting delta
for (var, lev) in zip([θsvd, θnmf], levels) 
    figure(figsize = (10, 8))
    anom = [v .- var[end] for v in var] 
    for (ii, i) in enumerate(inds)
        ax = subplot(3,3,ii, projection = ccrs.PlateCarree())
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
end

lab_pt = (360-51, 57)
Oũsvd = ptobserve(lab_pt, γ, spatialmodessvd, ũsvd) 
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

subplot(1,3,2)
Stimeseries = vec(Oũsvd.x[DimensionalData.Between(1920yr, 1980yr), At(:δ)]) ./ 0.5255
θtimeseries = vec(Oũsvd.x[DimensionalData.Between(1920yr, 1980yr), At(:θ)])

scatter(value.(Stimeseries), value.(θtimeseries), c = ustrip.(collect(Tᵤ[93:end])), cmap = "viridis", vmin = 1920, vmax = 1995)
[text(x, y, s) for (x,y,s) in zip(ustrip.(value.(Stimeseries)), ustrip.(value.(θtimeseries)), string.(convert.(Int64, ustrip.(Tᵤ[93:end]))))]
xlabel("Salinity anomaly from 1980 [g/kg]")
ylabel("θ anomaly from 1980 [K]")
title("SVD: T, inferred S from 51°W, 57°N")


subplot(1,3,3)
Stimeseries = vec(Oũnmf.x[DimensionalData.Between(1920yr, 1980yr), At(:δ)]) ./ 0.5255
θtimeseries = vec(Oũnmf.x[DimensionalData.Between(1920yr, 1980yr), At(:θ)])

scatter(value.(Stimeseries), value.(θtimeseries), c = ustrip.(collect(Tᵤ[93:end])), cmap = "viridis", vmin = 1920, vmax = 1995)
[text(x, y, s) for (x,y,s) in zip(ustrip.(value.(Stimeseries)), ustrip.(value.(θtimeseries)), string.(convert.(Int64, ustrip.(Tᵤ[93:end]))))]
xlabel("Salinity anomaly from 1980 [g/kg]")
ylabel("θ anomaly from 1980 [K]")
title("NMF: T, inferred S from 51°W, 57°N")
