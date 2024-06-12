
import Pkg
Pkg.activate("../")

using OPTinv, Unitful, DimensionalData, BLUEs, Statistics, UnitfulLinearAlgebra, LinearAlgebra, PyPlot, Revise, DrWatson, PyCall, TMI, CSV, DataFrames, JLD2, NaNMath, Dates, RollingFunctions, PaleoData, DateFormats
import Measurements.value as value
import OPTinv.Est

lab_pt = (360-51, 57)


ccrs = pyimport("cartopy.crs")
cm = pyimport("cmocean.cm")
cfeature = pyimport("cartopy.feature")
cu = pyimport("cartopy.util")
mpath = pyimport("matplotlib.path")
patches = pyimport("matplotlib.patches")

#which cores do you want?
corenums_full = [28, 26, 25, 22, 21, 20, 19, 10, 9,13,14]
core_list_full = Symbol.("MC".* string.() .* "A")
#missing  14
#corenums = [28, 26, 25, 22, 21, 20, 19, 10, 9,13,14]
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

@time yde, ỹdesvd, u₀de, ũsvd = solvesystem(y, u₀svd, Esvd, predictsvd);
@time yde, ỹdenmf, u₀de, ũnmf = solvesystem(y, u₀nmf, Enmf, predictnmf);
ỹdenmf

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
savefig(plotsdir("reconsol.png"))

figure(figsize = (8, 8))
for (i, m) in enumerate(ℳnmf.dims[2])
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
=#
llssvd = linearleastsquares(ustrip.(value.(vec(ũsvd.x[:, :, At(:δ)]))), ustrip.(value.(vec(ũsvd.x[:, :, At(:θ)]))))
println("SVD slope [K/permil] = " * string(llssvd[1]))
llsnmf = linearleastsquares(ustrip.(value.(vec(ũnmf.x[:, :, At(:δ)]))), ustrip.(value.(vec(ũnmf.x[:, :, At(:θ)]))))
println("NMF slope [K/permil] = " * string(llsnmf[1]))


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
anom_index = findall(x->1850yr < x < 1970yr, Tᵤ)
θsvdrs = cat([γreshape(θsvd[i], γ) for i in 1:length(θsvd)]..., dims = 3)
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
ax1 = gca()
clrs = ["red", "blue"]
labels = ["SVD", "NMF"] 
#for (cs, θ, lab) in zip(colors, [θsvdrs, θnmfrs], labels)
for (cs, θ, lab) in zip(clrs, [θsvdrs, θnmfrs], labels)
    θ = ustrip.(value.(θ))
    #θ .-= NaNMath.mean(vec(θ))
    θglobal = [NaNMath.mean(vec(θ[:, :, t])) for t in 1:length(Tᵤ)]
    θbox = boxmean(θ, γ.lat, γ.lon, 49, 89, 309, 21)
    #plot(DimArray(value.(θglobal .- mean(θglobal)), Ti(Tᵤ)), color = cs, label = lab * ": global")
    
    plot(DimArray(value.(θbox .- mean(θbox[anom_index])), Ti(Tᵤ)), color = cs, label = lab, linewidth = 5, zorder = 10000)#* ": 50°W-20°E, 50°N-90°N")
end

title("Mean Temperature Anomaly from 1850-1980yr: 50°W-20°E, 50°N-90°N") 

xlim([ustrip(minimum(T)), ustrip(T[end-1])])

lmrdataset = loadLMR("sst")
lmrsst = makeNaN(mean(lmrdataset["sst"][:, :, :, :], dims = 3)[:, :, 1, :])
lmrtime = lmrdataset["time"][:]; lmrlon = lmrdataset["lon"][:]; lmrlat = lmrdataset["lat"][:]
anomindex = findall(x->1850 < x < 1980, year.(lmrtime))
lmrsst = lmrsst .- mean(lmrsst[:, :, anomindex], dims = 3)[:, :, 1]


lmrtimeroll = rolling(mean, year.(lmrtime), 50)
lmrglobal = rolling(mean, [NaNMath.mean(lmrsst[:, :, i]) for i in 1:size(lmrsst)[3]], 50)

lmrbox = rolling(mean, boxmean(lmrsst, lmrlat, lmrlon, 50, 90, 360-50, 20), 50)
#plot(DimArray((lmrglobal .- mean(lmrglobal)) * K, Ti(lmrtime)), color = "black", label = "LMR: global", zorder = 0)
plot(DimArray((lmrbox  .- mean(lmrbox))* K, Ti(lmrtimeroll)), color = "purple", label = "LMR", zorder = 0, linewidth = 3)

hadisst = loadHadISST()
hadsst = makeNaN(hadisst["tos"][:, :, :])
hadlat = hadisst["latitude"][:]; hadlon = hadisst["longitude"][:]; hadtime = hadisst["time"][:]
hadroll(x) = rolling(mean, x, 50 * 12)
hadisst_gm = hadroll([NaNMath.mean(hadsst[:, :, i]) for i in 1:size(hadsst)[3]])
hadisst_bm = hadroll(boxmean(hadsst, hadlat, hadlon, 50, 90, -50, 20))
hadtimeroll = hadroll(yeardecimal.(hadtime))
#plot(DimArray(hadisst_gm .- hadisst_gm[begin], Ti(hadtimeroll)), color = "green", label = "HadISST: global")
plot(DimArray(hadisst_bm .- mean(hadisst_bm), Ti(hadtimeroll)), color = "green", label = "HadISST", linewidth = 3)

ylabel("Temperature Anomaly [K]")
xlabel("Time [yr CE]")

oc2k_binned, _, oc2k_ages = loadOcean2kBinned()
oc2k, _ = loadOcean2k()

atl_indices = findall(x->occursin("Atlantic", x), oc2k[!, "name"])
arc_indices = findall(x->occursin("Arctic", x), oc2k[!, "name"])
inregion = findall(x->x>50, oc2k[!, "lat"])
rel_records = intersect(vcat(atl_indices, arc_indices), inregion)

oc2kglobal = [NaNMath.mean(oc2k_binned[i, :]) for i in 1:size(oc2k_binned)[1]]
oc2kregion = [NaNMath.mean(oc2k_binned[i, rel_records]) for i in 1:size(oc2k_binned)[1]]
legend()
twin = ax1.twinx()
[twin.plot(oc2k_ages[:, i], oc2k_binned[:, i], color = "navy", alpha = 0.25, zorder = 1000, label = i == 1 ? "Ocean2k: Global" : "") for i in 1:size(oc2k_ages)[2]]
[twin.plot(oc2k_ages[:, i], oc2k_binned[:, i], color = "aqua", alpha = 0.25, zorder = 500, label = i == rel_records[1] ? "Ocean2k: Atl." : "") for i in rel_records]

binleft = 0:200:2000
[twin.plot([b, b + 200], [v, v], color = "navy", linewidth = 3, zorder = 0, label = b == binleft[1] ? "Ocean2k: Global Mean" : "") for (b, v) in zip(binleft, reverse(oc2kglobal))]
[twin.plot([b, b + 200], [v, v], color = "aqua", linewidth = 3, zorder = 0, label = b == binleft[1] ? "Ocean2k: Atlantic Mean" : "") for (b, v) in zip(binleft, reverse(oc2kregion))]
twin.set_ylabel("Standardized SST Anomaly [s.d. units]")

legend()
ax1.patch.set_visible(false)
ax1.set_zorder(twin.get_zorder()+1)
savefig(plotsdir("meants.png"))

res = 5
#plot 50 year means starting at 1500

inds =  findall(x->x == 1500.0yr, Tᵤ)[1]:res:length(Tᵤ)-res
#inds = inds[end-1:end]
@show length(inds)
lvls = [-1:0.1:1, -1:0.1:1, -1:0.1:1]
#stopped plotting delta
savenames = ["svd", "nmf", "lmrsst"]
proj = ccrs.Orthographic(central_longitude=-80+55,central_latitude = 35)
noproj = ccrs.PlateCarree(central_longitude = 0)

for (var, lev, sn) in zip([θsvd, θnmf, lmrsst], lvls, savenames)
#for (var, lev, sn) in zip([θnmf], lvls, savenames) 
    figure(figsize = (10, 10))
    
    for (ii, i) in enumerate(inds)
        ax = subplot(3,3,ii, projection = proj)
        ax_hdl = ax.plot(orthographic_axes(80,35,-80,30,5)...,
    color="black", linewidth=0.5,
            transform=noproj)
        tx_path = ax_hdl[1]._get_transformed_path()
        path_in_data_coords, _ = tx_path.get_transformed_path_and_affine()
        polygon1s = mpath.Path(path_in_data_coords.vertices)
        #vcode = [1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1]  #Path-code
        #polygon1v = mpath.Path(path_in_data_coords.vertices, vcode)
        ax.set_boundary(polygon1s)

        
        title(string(convert(Int64, ustrip(Tᵤ[i]))) * "-" *string(convert(Int64, ustrip(Tᵤ[i + res]))))
        ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "110m", edgecolor="k", facecolor="gray"))
        #ax.set_extent([-80, 30, 0, 90])
        #this is kind of a slow reshape but c'est la vie
        if sn != "lmrsst"
            anom = var #[v .- mean(var[anom_index]) for v in var] 
            plotme = γreshape(ustrip.(value.(sum(anom[i:i+res] ./ res))), γ)'
            #add a cyclic point 
            plotme, lon = cu.add_cyclic_point(plotme, coord = ustrip.(γ.lon))
            lat = γ.lat
        else
            
            tustart = ustrip(Tᵤ[i]); tuend = ustrip(Tᵤ[i + res])
            lmrindices = findall(x->tustart < x < tuend, year.(lmrtime))
            plotme = mean(lmrsst[:, :, lmrindices], dims = 3)[:, :, 1]'
            plotme, lon = cu.add_cyclic_point(plotme, coord = lmrlon)
            lat = lmrlat
        end
        
        cf = contourf(lon, lat, plotme, cmap = cm.balance, levels = lev, transform = noproj)
        c = contour(lon, lat, plotme, colors = "black", levels = lev, transform = noproj)
        ax.clabel(c)
        ax.gridlines(draw_labels = true)
    end
    tight_layout()
    savefig(plotsdir(sn * "surfacesol.png"))
end

 
lab_pt = (360-51, 57)
Oũsvd = ptobserve(lab_pt, γ, spatialmodessvd, ũsvd) 
Oũnmf = ptobserve(lab_pt, γ, spatialmodesnmf, ũnmf)

#en4time, en4temp, en4salinity = loadEN4("/home/brynn/Code/oceanFP/data/EN4/analyses/", lab_pt[1], lab_pt[2])
en4time, en4temp, en4salinity = loadEN4("/home/brynn/Code/oceanFP/data/EN4/analyses/", lab_pt[1], lab_pt[2])
#if you want to grab a box instead 
lab_lon = (360-60, 360-45)
lab_lat = (45, 64)
en4time, en4temp, en4salinity = loadEN4("/home/brynn/Code/oceanFP/data/EN4/analyses/", lab_lon, lab_lat)

df = CSV.read(datadir("TLS1998Fig7Digitized.csv"), DataFrame)
#figure(figsize = (16,4))
#subplot(1,3,1)
tls_anom_ind = findall(x->x == 70, df[!, "Year"])
tls_sal = df[!, "Salinity"] .- df[tls_anom_ind, "Salinity"]
tls_temp = df[!, "Temperature"] .- df[tls_anom_ind, "Temperature"]

en4timeann = unique(year.(en4time))
en4_anom_index = findall(x->x == 1970, en4timeann)
en4timeann = rolling(mean, en4timeann, 10)
en4_sal = monthlyannual(en4salinity, en4time)
en4_sal .-= en4_sal[en4_anom_index]
en4_sal = rolling(mean, en4_sal, 10)
en4_temp = monthlyannual(en4temp, en4time)
en4_temp .-= en4_temp[en4_anom_index]
en4_temp = rolling(mean, en4_temp, 10)

hadsst_lab = PaleoData.ptobserve(hadsst, hadlon, hadlat, lab_pt[1], lab_pt[2])
hadsst_lab = boxmean(hadsst, hadlat, hadlon, lab_lat[1], lab_lat[2], lab_lon[1], lab_lon[2])
hadtimeann = unique(year.(hadtime))
had_anom_ind = findall(x->x == 1970, hadtimeann)
hadsst_lab
hadsst_lab = monthlyannual(hadsst_lab, hadtime)
hadsst_lab .-= hadsst_lab[had_anom_ind]
hadtimeann = rolling(mean, hadtimeann, 10)
hadsst_lab = rolling(mean, hadsst_lab, 10)

figure(figsize = (8, 8))
subplot(2,1,1)
title("Labrador Sea T and S")
ou_anom_ind = findall(x->x == 1970yr, Tᵤ)
Ousvdθ = value.(ustrip.(Oũsvd.x[:, At(:θ)]))
Ousvdθ .-= Ousvdθ[ou_anom_index]
Ousvdδ = value.(ustrip.(Oũsvd.x[:, At(:δ)]))
Ousvdδ .-= Ousvdδ[ou_anom_ind] ./= 0.5255

Ounmfθ = value.(ustrip.(Oũnmf.x[:, At(:θ)]))
Ounmfθ .-= Ounmfθ[ou_anom_index]
Ounmfδ = value.(ustrip.(Oũnmf.x[:, At(:δ)]))
Ounmfδ .-= Ounmfδ[ou_anom_ind] ./= 0.5255

scatter(ustrip.(Tᵤ), Ousvdθ, color = "red", label = "SVD", s = 100, edgecolors = "black", zorder = 2)
scatter(ustrip.(Tᵤ), Ounmfθ, color = "blue", label = "NMF", s = 100, edgecolors = "black", zorder = 2)
plot(1900 .+ df[!, "Year"], tls_temp, color = "black", label = "TLS1998", zorder = 1)
plot(en4timeann, en4_temp, color = "gray", label = "EN4", zorder = 1)
plot(hadtimeann, hadsst_lab, color = "teal", label = "HadiSST", zorder = 1)
legend()
xlim(1850, 1970)
ylim(-1.5,1.5)
ylabel("Temp. Anom. from 1970 [K]")

subplot(2,1,2)
plot(1900 .+ df[!, "Year"], tls_sal, color = "black", label = "TLS1998")
plot(en4timeann, en4_sal, color = "gray", label = "EN4")

scatter(ustrip.(Tᵤ), Ousvdδ , color = "red", label = "SVD", s = 100, edgecolors = "black", zorder = 1)
scatter(ustrip.(Tᵤ), Ounmfδ, color = "blue", label = "NMF", s = 100, edgecolors = "black", zorder = 1)
xlim(1850, 1970)
ylim(-0.3, 0.3)
ylabel("Sal. Anom. from 1970 [g/kg]")
xlabel("Time [years CE]") 
tight_layout()
#plot(tls_sal, tls_temp, c = df[!, "Year"], cmap = "viridis", marker = "s", edgecolors = "purple", s = 100)
#colorbar(s)
#scatter(en4_sal, en4_temp, c = en4timeann, cmap = "viridis", vmin = 1920, vmax = 1995, marker = "*", edgecolors = "black", s = 100)
#plot(en4_sal, en4_temp)
#xlabel("Salinity Anomaly from 1970")
#ylabel("Temperature Anomaly from 1970")

#[text(x,y,s) for (x,y,s) in zip(df[!, "Salinity"], df[!, "Temperature"], string.(df[!, "Year"])), vmin = 20, vmax = 95]
plot(df[!, "Salinity"], df[!, "Temperature"], color = "gray") 
xlabel("Salinity")
ylabel("Temperature") 
title("Digitized TLS1998 Figure 7")

subplot(1,3,2)
Stimeseries = vec(Oũsvd.x[DimensionalData.Between(1920yr, 1970yr), At(:δ)]) ./ 0.5255
θtimeseries = vec(Oũsvd.x[DimensionalData.Between(1920yr, 1970yr), At(:θ)])
t = collect(1920:10:1970)
scatter(value.(Stimeseries), value.(θtimeseries), c = t, cmap = "viridis", vmin = 1920, vmax = 1995)
[text(x, y, s) for (x,y,s) in zip(ustrip.(value.(Stimeseries)), ustrip.(value.(θtimeseries)), string.(t))]
xlabel("Salinity anomaly from 1980 [g/kg]")
ylabel("θ anomaly from 1980 [K]")
title("SVD: T, inferred S from 51°W, 57°N")


subplot(1,3,3)
Stimeseries = vec(Oũnmf.x[DimensionalData.Between(1920yr, 1970yr), At(:δ)]) ./ 0.5255
θtimeseries = vec(Oũnmf.x[DimensionalData.Between(1920yr, 1970yr), At(:θ)])

scatter(value.(Stimeseries), value.(θtimeseries), c = t, cmap = "viridis", vmin = 1920, vmax = 1995)
[text(x, y, s) for (x,y,s) in zip(ustrip.(value.(Stimeseries)), ustrip.(value.(θtimeseries)), string.(t))]
xlabel("Salinity anomaly from 1980 [g/kg]")
ylabel("θ anomaly from 1980 [K]")
title("NMF: T, inferred S from 51°W, 57°N")
savefig(plotsdir("labcomp.png"))

#=
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
    ptsvd = ptobserve((lon, lat), γ, spatialmodessvd, ũsvd)
    sst = convertx(data[!, name * "_SST"]) 
    age = convertx(data[!, name * "_Age"])
    subplot(4,2, ii)
    plot(1950 .- age, sst .- NaNMath.mean(sst), color = "black")
    plot(value.(ptnmf.x[:, At(:θ)]), color = "blue")
    plot(value.(ptsvd.x[:, At(:θ)]), color = "red")
    title(name)
    xlim(ustrip(minimum(T)),1970)
    ylim(-2, 2)
    xlabel("Time [years CE]")
    ylabel("Temperature anomaly [K]")
end
tight_layout()
savefig(plotsdir("ocean2k.png"))

f = figure()
ax = f.add_subplot(projection = ccrs.PlateCarree())
lats = meta[intersect(inregion, atl_indices), "lat"]
lons = meta[intersect(inregion, atl_indices), "lon"]
names = meta[intersect(inregion, atl_indices), "name"]
ax.coastlines(color = "gray")
ax.set_extent([-30, 30, 50, 70])
scatter(lons, lats)
[text(lon, lat, t[13:end]) for (lon, lat, t) in zip(lons, lats, names)]
savefig(plotsdir("ocean2klocs.png"))
=#

WHtempnmf = ptobserve((360-45, 57), γ, spatialmodesnmf, ũnmf)
WHtempsvd = ptobserve((360-45, 57), γ, spatialmodessvd, ũsvd)
EGtempnmf = ptobserve((360-14, 75), γ, spatialmodesnmf, ũnmf)
EGtempsvd = ptobserve((360-14, 75), γ, spatialmodessvd, ũsvd)

figure();plot(value.(WHtempnmf.x[:, At(:θ)]), color = "red", label = "S. Greenland Anom.")
plot(value.(WHtempsvd.x[:, At(:θ)]), color = "red", linestyle = "dashed")
xlim(ustrip(minimum(T)), 1970)

plot(value.(EGtempnmf.x[:, At(:θ)]), color = "blue", label = "E. Greenland Anom.")
plot(value.(EGtempsvd.x[:, At(:θ)]), color = "blue", linestyle = "dashed")
legend()
ylim(-1.5, 1.5)
twin = gca().twinx()
thornalley = loadThornalley2018()

for k in keys(thornalley)
    age = thornalley[k][!, "age [CE]"]
    data = thornalley[k][!, "smooth"]
    twin.plot(age, data .- NaNMath.mean(data), color = "black", alpha = 0.5, zorder = 0)
end
ylabel(L"\bar{SS}" * " [mm]")
savefig(plotsdir("EGWG.png"))

