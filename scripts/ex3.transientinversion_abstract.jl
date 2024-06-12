
import Pkg
Pkg.activate("../")

using OPTinv, Unitful, DimensionalData, BLUEs, Statistics, UnitfulLinearAlgebra, LinearAlgebra, PyPlot, Revise, DrWatson, PyCall, TMI, CSV, DataFrames, JLD2, NaNMath, Dates, RollingFunctions, PaleoData, DateFormats, Measurements
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

shallow = Symbol.("MC" .* string.([28, 26, 25]) .* "A")
intermediate = Symbol.("MC" .* string.([22, 21, 20, 19]) .* "A")
deep = Symbol.("MC" .* string.([10, 9,13,14]) .* "A")


@time allc = invert(inversion(allcores, "svd", "all", "red"))
@time oldc = invert(inversion(oldcores, "svd", "old", "blue"))
                        
solutions = [oldc, allc]
suffix = "allvold"

# @time shallowc = invert(inversion(shallow, "svd", "shallow", "purple"))
# @time intermediatec = invert(inversion(intermediate, "svd", "intermediate", "green"))
# @time deepc = invert(inversion(deep, "svd", "deep", "orange"))

#solutions = [shallowc, intermediatec, deepc]
#suffix = "depth"

# ============ FIGURE 1:SOLUTIONS RECONSTRUCTED AT CORE SITES ===================== #
locs = core_locations() 
min_age = [minimum(s.y.dims[1]) for s in solutions]
max_age = [maximum(s.y.dims[1]) for s in solutions]
fig = figure(figsize = (8, 8))
xl = [minimum(min_age), maximum(max_age)]
for (i, core) in enumerate(allcores) 
    subplot(4,3, i)
    #we want to plot y from the solution that has the longest record 
    @show core_in_sol = [core ∈ s.y.dims[2] for s in solutions]
    ma_mod = copy(min_age)
    ma_mod[(!).(core_in_sol)] .= 1980yr
    @show ma_mod
    @show plotfrom = findmin(ma_mod)[2]
    for (j, s) in enumerate(solutions)
        if j == plotfrom
            plot(s.y.x[:, At(core)], label = "y", color = "black", zorder = j)
        end
        if core in s.y.dims[2]
            plot(s.ỹ.x[:, At(core)], label = "ỹ: " * s.name, color = s.color, zorder = j)
        end
    end
    
    ylim(-0.25, 0.25)
    xlim(ustrip.(xl))


    if i ∈ [9,10,11]
        xticks(1000:250:1750, fontsize = 12)
        xlabel("Time [years CE]", fontsize = 15)
    else
        xticks(1000:250:1750, ["", "", "", ""])
        xlabel("")
    end
    
    if i ∈ [1,4,7,10]
        yticks(-0.2:0.2:0.2, fontsize = 12) 
        ylabel(L"\delta^{18}\mathrm{O}_\mathrm{calcite}" * " [‰]", fontsize = 15)
    else
        ylabel("")
        yticks(-0.2:0.2:0.2, ["", "", ""])
    end
        
    title(string(core) * ", " * string(locs[core][3]) * "m", fontsize = 15)
    grid()
    gca().invert_yaxis()
end
tight_layout()
savefig(plotsdir("reconsol" * suffix * ".png"))

# ================FIGURE 2: MODE MAGNITUDES============== #
for s in [:θ, :δ]
figure(figsize = (8,8))
    for (i, sol) in enumerate(solutions)
        subplot(length(solutions),1, i)
        modes = sol.ũ.dims[2]
        for m in 1:3
            #color = cmap.get_cmap("nipy_spectral")(1/(maximum(modes) - minimum(modes)) * m)
            color = ["purple", "green", "darkorange"][m]
            ls = ["solid", "dashed", "dotted"][m]
            tr = (minimum(sol.y.dims[1]), maximum(sol.y.dims[1]))
            ts = value.(sol.ũ.x[:, At(m), At(s)])
            plot((ts .- mean(ts)) ./ std(ts), label = string(m), color = color, zorder = length(modes) - m, linewidth = 4, linestyle = ls)
            
        end
        if i == 2

            xticks(600:200:2000, fontsize = 12)
            xlabel("Time [years CE]", fontsize = 15)
        else
            xticks(600:200:2000, fill("", length(600:200:2000)))
            xlabel("")
        end
        
        yticks(-4:2:6, fontsize = 12) 
        ylabel("Mode Magnitude Anomaly [σ]", fontsize = 15)
        #legend(loc = "center left")
        text(x = 485, y = 5, s = ["A", "B"][i], fontsize = 30, fontweight = "bold")
        grid()
        
        xlim(475, 2000)
        ylim(-4, 6)
    
    end
    tight_layout()
    savefig(plotsdir("modemags"*string(s)* suffix* ".png"))

end

# ================= FIGURE 3: REGION MEAN ====================== #
regionmeanindices = γbox(oldc.γ, 49, 89, 309, 21)
figure()
grid()
for sol in solutions
    sol_anom_ind = findall(x->1850yr<x<1970yr, Array(sol.ũ.dims[1]))
    θbox = [mean(sol.θ[i][regionmeanindices]) for i in 1:length(sol.θ)]
    θbox .-= mean(θbox[sol_anom_ind])
    inds = findall(x->x∈sol.y.dims[1], Array(sol.ũ.dims[1]))
    
    @show min = findmin(θbox)
    @show sol.ũ.dims[1][min[2]]
    pre1900inds = intersect(inds, findall(x->1000yr<x<1900yr, Array(sol.ũ.dims[1])))
    @show max = findmax(θbox[pre1900inds])
    @show sol.ũ.dims[1][pre1900inds][max[2]]

    
    @show (max[1] - min[1]) / (sol.ũ.dims[1][min[2]] - sol.ũ.dims[1][pre1900inds][max[2]]) *1000

    println("LLS")
    maxind = findall(x->x == sol.ũ.dims[1][pre1900inds][max[2]], Array(sol.ũ.dims[1]))[1]
    x = ustrip.(sol.ũ.dims[1][maxind:min[2]])
    
    C = Diagonal(Measurements.uncertainty.(θbox[maxind:min[2]]).^2)
    y = ustrip.(value.(θbox[maxind:min[2]]))
    lls, C = linearleastsquares(x, y; C)
    @show lls[1] * 1000 
    @show sqrt.(diag(C))[1] .* 1000
    plot(DimArray(value.(θbox[inds]), sol.ũ.dims[1][inds]), color = sol.color, label = sol.name, linewidth = 3, zorder = 10000)
    println()
end


#title("Mean Temperature Anomaly from 1850-1970yr: 50°W-20°E, 50°N-90°N") 
lmrdataset = loadLMR("sst")
lmrtime = lmrdataset["time"][:]; lmrlon = lmrdataset["lon"][:]; lmrlat = lmrdataset["lat"][:]
lmrsst = makeNaN(mean(lmrdataset["sst"][:, :, :, :], dims = 3)[:, :, 1, :])
anomindex = findall(x->1850 < x < 1970, year.(lmrtime))
lmrsst = lmrsst .- mean(lmrsst[:, :, anomindex], dims = 3)[:, :, 1]

lmrtimeroll = rolling(mean, year.(lmrtime), 50)
lmrglobal = rolling(mean, [NaNMath.mean(lmrsst[:, :, i]) for i in 1:size(lmrsst)[3]], 50)

lmrbox = rolling(mean, boxmean(lmrsst, lmrlat, lmrlon, 50, 90, 360-50, 20), 50)
#plot(DimArray((lmrglobal .- mean(lmrglobal)) * K, Ti(lmrtime)), color = "black", label = "LMR: global", zorder = 0)
plot(DimArray((lmrbox  .- mean(lmrbox))* K, Ti(lmrtimeroll)), color = "darkgray", label = "LMR", zorder = 0, linewidth = 3, linestyle = "dashdot")

hadisst = loadHadISST()
hadsst = makeNaN(hadisst["tos"][:, :, :])
hadlat = hadisst["latitude"][:]; hadlon = hadisst["longitude"][:]; hadtime = hadisst["time"][:]
hadroll(x) = rolling(mean, x, 50 * 12)
hadisst_gm = hadroll([NaNMath.mean(hadsst[:, :, i]) for i in 1:size(hadsst)[3]])
hadisst_bm = hadroll(boxmean(hadsst, hadlat, hadlon, 50, 90, -50, 20))
hadtimeroll = hadroll(yeardecimal.(hadtime))
#plot(DimArray(hadisst_gm .- hadisst_gm[begin], Ti(hadtimeroll)), color = "green", label = "HadISST: global")
plot(DimArray(hadisst_bm .- mean(hadisst_bm), Ti(hadtimeroll)), color = "black", label = "HadISST", linewidth = 3, linestyle = "dashed")

ylabel("Temperature Anomaly [K]", fontsize = 15)
xlabel("Time [yr CE]", fontsize = 15)

oc2k_binned, oc2k_binnedv, oc2k_ages = loadOcean2kBinned()
oc2k, oc2kdata = loadOcean2k()

atl_indices = findall(x->occursin("Atlantic", x), oc2k[!, "name"])
arc_indices = findall(x->occursin("Arctic", x), oc2k[!, "name"])
inregion = findall(x->x>50, oc2k[!, "lat"])
rel_records = intersect(vcat(atl_indices, arc_indices), inregion)

oc2kglobal = [NaNMath.mean(oc2k_binned[i, :]) for i in 1:size(oc2k_binned)[1]]
oc2kregion = [NaNMath.mean(oc2k_binned[i, rel_records]) for i in 1:size(oc2k_binned)[1]]


binleft = 0:200:2000
bincenter = 100:200:2000
binright = 200:200:2000
size(oc2k)
oc2kanom = Matrix{Float64}(undef, size(oc2k)[1], length(bincenter))
oc2kanom .= NaN
for i in 1:57
    sstname = names(oc2kdata)[i*2-1]
    if occursin("Atlantic", sstname) 
        meta_index = findall(x->occursin(sstname[begin:begin+11], x), oc2k.name)[1]
        lat = oc2k[meta_index, "lat"]
            t = oc2kdata[!, sstname]
            t[t .== "NaN"] .= NaN
            y = oc2kdata[!, names(oc2kdata)[i*2]]
            y[y .== "NaN"] .= NaN
            for (j, (bl, br)) in enumerate(zip(binleft, binright))
                ind = findall(x->bl<x<br, t)
                if length(ind) > 0
                    oc2kanom[i, j] = mean(y[ind])
                else
                    oc2kanom[i, j] = NaN
                end
            end
    end
    
end

oc2kanom .-= oc2kanom[:, end]
#[plot(bincenter, oc2kanom[i, :], color = "gray") for i in 1:57]
#plot(bincenter, [NaNMath.mean(oc2kanom[:, i]) for i in 1:10], color = "black", label = "Ocean2k: Atlantic Mean", ".-")
slope_est = (0.38, 0.54)
t = Array(ustrip.(oldc.ũ.dims[1]))
sol_anom_ind = findall(x->1850<x<1970, t)
t_oc2k = 1500yr:1yr:1700yr
y1 = [(1970 - t_) * slope_est[1] * 0.001 for t_ in ustrip.(t_oc2k)]
y2 = [(1970 - t_) * slope_est[2] * 0.001 for t_ in ustrip.(t_oc2k)]
#fill_between(x = t, y1 = y1 .- mean(y1[sol_anom_ind]), y2 = y2 .- mean(y2[sol_anom_ind]), color = "black", alpha = 0.5, label = "Global Ocean2k SST Trend")
plot(ustrip.(t_oc2k), y1 .+ 0.05, color = "black", label = "Global Oc2k SST Trend")
text(x = 1709, y = 0.135, s = "0.38K/kyr")
plot(ustrip.(t_oc2k), y2 .+ 0.05, color = "black")
text(x = 1709, y = 0.18, s = "0.54K/kyr")
legend()

Tm1 = ustrip.(minimum([sol.y.dims[1][1] for sol in solutions]))
Tm2 = ustrip.(maximum([sol.y.dims[1][end-1] for sol in solutions]))
xlim(Tm1, Tm2)
xticks(800:200:1800, fontsize = 12) 
ylim(-0.5, 0.5)
yticks(-0.5:0.25:0.5, fontsize = 12)
tight_layout()
savefig(plotsdir("meants" * suffix * ".png"))

# ================= FIGURE 4: MAPS ============= #
inds = collect(1100yr:200yr:2000yr)
stepsize = 70 
lev = -1.5:0.2:1.5

proj = ccrs.Orthographic(central_longitude=-80+55,central_latitude = 35)
noproj = ccrs.PlateCarree(central_longitude = 0)
fig, axs = subplots(nrows = length(solutions) + 1, ncols = length(inds),
                    subplot_kw = Dict("projection"=>proj), figsize = (15, 8))
#(80,35,-80,30,5)
lon_box = [-80, 30]
lat_box = [35, 80]
for (i, sol) in enumerate([solutions..., lmrsst])
    Tu = sol isa solution ? Array(sol.ũ.dims[1]) : year.(lmrtime)yr
    Ty = sol isa solution ? Array(sol.y.dims[1]) : year.(lmrtime)yr 
    anom_index = findall(x->1850yr<x<1970yr, Tu)

    #remove 1850-1970anomaly from value we are plotting 
    if sol isa solution
        θ_ = [ustrip.(value.(x)) for x in sol.θ]
        θ = [x .- (sum(θ_[anom_index]) ./ length(anom_index)) for x in θ_]
    else
        θ = sol
        θ .-= (sum(θ[:, :, anom_index], dims = 3) ./ length(anom_index))
    end
            
    for (ii, y) in enumerate(inds)
        
        start_index = findall(x->x == y, Tu)[1]
        stop_index = findall(x->x == y + stepsize*yr, Tu)[1]
    
        if sol isa solution && y ∈ Ty
            plotme = γreshape(sum(θ[start_index:stop_index] ./ (stepsize/10)), sol.γ)' #100yr at 10yr res 
        elseif sol isa Array
            plotme = sum(θ[:, :, start_index:stop_index], dims = 3)[:, :, 1] ./ stepsize #100yr at 1yr res
            plotme = plotme'
        end
        if (sol isa solution && y ∈ Ty) || (sol isa Array)
            ax =  axs[(ii-1) * (length(solutions)+1)+ i]
            ax_hdl = ax.plot(orthographic_axes(lat_box..., lon_box..., 5)...,
                             color="black", linewidth=0.5,
                             transform=noproj)
            tx_path = ax_hdl[1]._get_transformed_path()
            path_in_data_coords, _ = tx_path.get_transformed_path_and_affine()
            polygon1s = mpath.Path(path_in_data_coords.vertices)
            ax.set_boundary(polygon1s)
            ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "110m", edgecolor="k", facecolor="gray"))
            
            ax.set_title(string(convert(Int64, ustrip(Tu[start_index]))) * "-" *string(convert(Int64, ustrip(Tu[stop_index]))))
            
            #add a cyclic point
            lon = sol isa solution ? ustrip.(sol.γ.lon) : lmrlon 
            #plotme, lon_ = cu.add_cyclic_point(plotme, coord = lon)
            lat = sol isa solution ? ustrip.(sol.γ.lat) : lmrlat
            lat_box_= [lat_box[1] - 20, lat_box[2] + 20]
            lon_box_ = [360-80 - 20, 30 + 20]
            plotme = geospatialsubset(plotme', lat, lon, lat_box_, lon_box_)
            plotme = plotme'
            cf = ax.contourf(lon, lat, plotme, cmap = cm.balance, levels = lev, transform = noproj)
            c = ax.contour(lon, lat, plotme, colors = "black", levels = lev, transform = noproj)
            ax.clabel(c)
            gl = ax.gridlines(draw_labels = true)
            gl.top_labels = false
            gl.right_labels = false
            gl.bottom_labels = false
        end
    end
    tight_layout()
    savefig(plotsdir("surfacesol" * suffix * ".png"))
end
# ===================== FIGURE 5alt: Labrador T and S  ================ #

df = CSV.read(datadir("TLS1998Fig7Digitized.csv"), DataFrame)

tls_anom_ind = findall(x->x == 70, df[!, "Year"])
tls_sal = df[!, "Salinity"] .- df[tls_anom_ind, "Salinity"]
tls_temp = df[!, "Temperature"] .- df[tls_anom_ind, "Temperature"]
tls_time = 1900 .+ df[!, "Year"]
lab_lon = (360-60, 360-45)
lab_lat = (55, 64)
labθ = [[mean(sol.θ[i][γbox(sol.γ, lab_lat[1], lab_lat[2], lab_lon[1], lab_lon[2])]) for i in 1:length(sol.θ)] for sol in solutions]
labS =[[mean(sol.δ[i][γbox(sol.γ, lab_lat[1], lab_lat[2], lab_lon[1], lab_lon[2])]) ./ 0.5255 for i in 1:length(sol.θ)] for sol in solutions]

tls_btime = 1920:10:1970
tls_bsal = Vector{typeof(tls_sal[1])}(undef, length(tls_btime))
tls_btemp = Vector{typeof(tls_temp[1])}(undef, length(tls_btime))
for (i,t) in enumerate(tls_btime)
    ind = findall(x->t<x<t+10, tls_time)
    tls_bsal[i] = mean(tls_sal[ind]) 
    tls_btemp[i] = mean(tls_temp[ind])
end

figure(figsize = (15, 5))
#[arrow(x = xi, y = yi, dx = dxi, dy = dyi, head_length = 0.01, head_width = 0.001,linewidth = 0.01, facecolor = "k", ec = "k") for (xi, yi, dxi, dyi) in zip(tls_bsal[begin:end-1], tls_btemp[begin:end-1], diff(tls_bsal), diff(tls_btemp))]

for (i, (sol, Ouθ, OuS, cmap)) in enumerate(zip(solutions, labθ, labS, ["Blues", "Reds"]))
    subplot(1, length(solutions), i) 
    s = scatter(x = tls_bsal .- tls_bsal[end], y = tls_btemp .- tls_btemp[end], c = tls_btime, cmap = "Greys", vmin = 1920, vmax = 1970, edgecolors = "black", s  = 150)
    plot(tls_bsal .- tls_bsal[end], tls_btemp .- tls_btemp[end], color = "black") 
    colorbar(s) 
    Tᵤ = Array(sol.ũ.dims[1])
    ou_anom_ind = findall(x->x == 1970yr, Tᵤ)
    ind = findall(x->x ∈ tls_btime .* yr, Tᵤ)
    Ouθ_ = value.(ustrip.(Ouθ))
    Ouθ_ .-= Ouθ_[ou_anom_ind]
    OuS_ = value.(ustrip.(OuS))
    OuS_ .-= OuS_[ou_anom_ind]
    
    s = scatter(OuS_[ind], Ouθ_[ind], c = ustrip.(Tᵤ)[ind], cmap = cmap, label = sol.name, s = 100, edgecolors = "white", vmin = 1920, vmax = 1970)
    plot(OuS_[ind], Ouθ_[ind], color = sol.color)
    colorbar(s)
    xlabel("Sal. Anom. [psu]")
    ylabel("Temp. Anom. [K]")
    xlim(-0.055, 0.055)
    ylim(-0.55, 0.55)
end
tight_layout()

savefig(plotsdir("labcomp_tands_" * suffix * ".png")) 


# ===================== FIGURE 5: Labrador T and S  ================ #

lab_pt = (360-51, 57)
#if you want to grab a box instead 
lab_lon = (360-60, 360-45)
lab_lat = (55, 64)

#Vector of DimEstimate
#pt observe
#Oũ = [OPTinv.ptobserve(lab_pt, sol.γ, sol.spatialmodes, sol.ũ) for sol in solutions]
#lab sea observe
labθ = [[mean(sol.θ[i][γbox(sol.γ, lab_lat[1], lab_lat[2], lab_lon[1], lab_lon[2])]) for i in 1:length(sol.θ)] for sol in solutions]
labS =[[mean(sol.δ[i][γbox(sol.γ, lab_lat[1], lab_lat[2], lab_lon[1], lab_lon[2])]) ./ 0.5255 for i in 1:length(sol.θ)] for sol in solutions]

en4time, en4temp, en4salinity = loadEN4("/home/brynn/Code/oceanFP/data/EN4/analyses/", lab_pt[1], lab_pt[2])
en4time, en4temp, en4salinity = loadEN4("/home/brynn/Code/oceanFP/data/EN4/analyses/", lab_lon, lab_lat)
mldtemp, mldsalinity, lon_rs, lat_rs, recs, recstime = loadEN4MLD();

df = CSV.read(datadir("TLS1998Fig7Digitized.csv"), DataFrame)
df = df[sortperm(df[!, "Year"]), :]

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


en4mld_temp = boxmean(mldtemp, lat_rs, lon_rs, lab_lat[1], lab_lat[2], -(360 - lab_lon[1]), -(360 - lab_lon[2]))
en4mld_time = 1900:2022
en4mld_temp = rolling(mean, en4mld_temp, 10)
en4mld_time = rolling(mean, en4mld_time, 10)
en4mld_anom = findall(x->x == 1970.5, en4mld_time)
en4mld_temp .-= en4mld_temp[en4mld_anom]

en4mld_sal =  boxmean(mldsalinity, lat_rs, lon_rs, lab_lat[1], lab_lat[2], -(360 - lab_lon[1]), -(360 - lab_lon[2]))
en4mld_sal = rolling(mean, en4mld_sal, 10) 
en4mld_sal .-= en4mld_sal[en4mld_anom]


hadsst_lab = PaleoData.ptobserve(hadsst, hadlon, hadlat, lab_pt[1], lab_pt[2])
hadsst_lab = boxmean(hadsst, hadlat, hadlon, lab_lat[1], lab_lat[2], lab_lon[1], lab_lon[2])
hadtimeann = unique(year.(hadtime))
had_anom_ind = findall(x->x == 1970, hadtimeann)
hadsst_lab = monthlyannual(hadsst_lab, hadtime, 1:3)
hadsst_lab .-= hadsst_lab[had_anom_ind]

#hadtimeann = rolling(mean, hadtimeann, 10)
#hadsst_lab = rolling(mean, hadsst_lab, 10)

figure(figsize = (8, 8))
subplot(2,1,1)
text(x = 1901, y = 0.75, s = "A", fontsize = 30, weight = "bold")
#title("Labrador Sea T and S")
grid()
for (sol, Ou) in zip(solutions, labθ)
    Tᵤ = Array(sol.ũ.dims[1])
    ou_anom_ind = findall(x->x == 1970yr, Tᵤ)
    ind = findall(x->x ∈ Array(sol.y.dims[1]), Tᵤ)
    Ouθ = value.(ustrip.(Ou))
    Ouθ .-= Ouθ[ou_anom_ind]
    scatter(ustrip.(Tᵤ)[ind], Ouθ[ind], color = sol.color, label = sol.name, s = 100, edgecolors = "black", zorder = 2)
end
    
plot(1900 .+ df[!, "Year"], tls_temp, color = "black", label = "TLS1998", zorder = 1, linestyle = "solid", linewidth = 3)
#plot(en4timeann, en4_temp, color = "black", label = "EN4 SST", zorder = 1, linestyle = "dotted")
#plot(en4mld_time, en4mld_temp, color = "black", label = "EN4 MLD", zorder = 1, linestyle = "dotted", linewidth = 3)
plot(hadtimeann, hadsst_lab, color = "black", label = "HadiSST", zorder = 1, linestyle = "dashed")
legend()
xlim(1900, 1970)
ylim(-1,1)
xticks(1900:10:1970, fontsize = 12)
yticks(-1:0.5:1, fontsize = 12)

ylabel("Temperature Anomaly [K]", fontsize = 15)



subplot(2,1,2)
text(x = 1901, y = 0.045, s = "B", fontsize = 30, weight = "bold")
for (sol, Ou) in zip(solutions, labS)
    Tᵤ = Array(sol.ũ.dims[1])
    ou_anom_ind = findall(x->x == 1970yr, Tᵤ)
    ind = findall(x->x ∈ Array(sol.y.dims[1]), Tᵤ)
    Ouδ = value.(ustrip.(Ou))
    Ouδ .-= Ouδ[ou_anom_ind]
    scatter(ustrip.(Tᵤ)[ind], Ouδ[ind], color = sol.color, label = sol.name, s = 100, edgecolors = "black", zorder = 2)
end

plot(1900 .+ df[!, "Year"], tls_sal, color = "black", label = "TLS1998", linestyle = "solid", linewidth = 3)
#plot(en4timeann, en4_sal, color = "black", label = "EN4 SSS", linestyle = "dotted")
#plot(en4mld_time, en4mld_sal, color = "black", label = "EN4 MLD", linestyle = "dotted", linewidth = 3)
xlim(1900, 1970)
grid()
ylim(-0.06, 0.06)
ylabel("Salinity Anomaly [g/kg]", fontsize = 15)
xlabel("Time [years CE]", fontsize = 15)
xticks(1900:10:1970, fontsize = 12)
yticks(-0.06:0.03:0.06, fontsize = 12)
tight_layout()
savefig(plotsdir("labcomp" * suffix * ".png")) 

## ================== FIGURE 6: E. Gr. and S. Gr. Anom. Mag.  ================= #
nordicind = γbox(oldc.γ, 65, 80,360-20, 20)
spg = γbox(oldc.γ, 50,65, 300,0)
spgW= γbox(oldc.γ, 50,65, 300,330)
spgE= γbox(oldc.γ, 50,65, 330,0)

inds = NamedTuple{(:Nordic, :SPG, :SPGE, :SPGW)}([nordicind, spg, spgE, spgW])

# WHtemp = [OPTinv.ptobserve((360-45, 57), sol.γ, sol.spatialmodes, sol.ũ) for sol in solutions]
# EGtemp = [OPTinv.ptobserve((360-14, 75), sol.γ, sol.spatialmodes, sol.ũ) for sol in solutions]
# SItemp = [OPTinv.ptobserve((360-20, 58), sol.γ, sol.spatialmodes, sol.ũ) for sol in solutions]

fig = figure(figsize = (10,8));
fig.subplots_adjust(right=0.75)
for (i, sol) in enumerate(solutions)
    
    ax = subplot(length(solutions), 1, i)
    grid()
    T = sol.y.dims[1]
    for (k, ls, lw, col) in zip(keys(inds), ["solid", "solid", "dotted", "dashed"], [4,4,3,3], ["purple", "darkgreen", "darkgreen", "darkgreen"])
        Tu = sol.ũ.dims[1]
        temp = DimArray([mean(sol.θ[i][inds[k]]) for i in 1:length(sol.θ)], Ti(Array(Tu)))
        anom = mean(ustrip.(value.(temp[DimensionalData.Between(1850yr, 1970yr)])))
        y1 = ustrip.(value.(temp[DimensionalData.Between(T[1], T[end])])) .- anom 
        plot(ustrip.(T), y1, color = col, label = k, linestyle = ls, linewidth = lw)
        @show col
        @show min = findmin(temp[DimensionalData.Between(T[1], T[end])])
        @show y1.dims[1][min[2]]
        @show max = findmax(temp[DimensionalData.Between(T[1], T[end])])
        @show y1.dims[1][max[2]]
        @show max[1] - min[1] 
        @show (max[1] - min[1]) / (y1.dims[1][min[2]] - y1.dims[1][max[2]]) * 1000
        println()
    end
    #title(sol.name)
    xlim(Tm1, Tm2)
    legend(loc = "upper right")
    ylabel("Temperature Anomaly [K]",fontsize = 15)
    if i == 2 xlabel("Time [years CE]", fontsize = 15) end
 
    yticks(-0.5:0.25:0.75, fontsize = 12)
    ylim(-0.5, 0.8)
    xticks(800:200:1800, fontsize = 12)
    text(x = 810, y = -0.48, s = ["A", "B"][i], fontsize = 30, weight = "bold")
    
    twin1 = gca().twinx()
    dahljensen = CSV.read(datadir("DigitizedDahlJensen.csv"), DataFrame)
    djyears = dahljensen[!, "Time [years CE]"]
    anom_index = findall(x->1850 < x < 1970, djyears)
    djtemp = dahljensen[!, "Temperature [degC]"] #.- mean(dahljensen[anom_index, "Temperature [degC]"])  
    twin1.plot(djyears, djtemp, color = "blue", linewidth = 1, alpha = 0.5)
    twin1.set_ylabel("Borehole Temp. [°C]", fontsize = 15, color = "blue")
    twin1.set_yticks(labels = -32.5:0.5:-30, ticks = -32.5:0.5:-30, fontsize = 12, color = "blue")
    
    
    twin2 = ax.twinx()
    twin2.spines.right.set_position(("axes", 1.15))
    steinhilber = loadSteinhilber2009()
    st_time = 1950 .- steinhilber[!, "YearBP"]
    st = Measurements.measurement.(steinhilber[!, "dTSI"], steinhilber[!, "dTSI_sigma"])
    twin2.plot(st_time, steinhilber[!, "dTSI"], color = "grey")
    twin2.fill_between(x = st_time, y1 = steinhilber[!, "dTSI"] .-  steinhilber[!, "dTSI_sigma"], y2 = steinhilber[!, "dTSI"] .+  steinhilber[!, "dTSI_sigma"], color = "grey", alpha = 0.2)
    twin2.set_ylabel("Total Solar Insolation [Wm⁻²]", fontsize = 15, color = "grey")
    twin2.set_yticks(labels = -1.5:0.5:1.5, ticks = -1.5:0.5:1.5, fontsize = 12, color = "grey")

    twin3 = ax.twinx()
    twin3.spines.right.set_position(("axes", 1.3))
    gao = loadGao2008()
    twin3.plot(gao[!, "Year"], gao[!, "Global"], color = "red", alpha = 0.5, linewidth = 1)
    twin3.set_ylabel("Global strat. sulf. aer. inj. [Tg]", fontsize = 15, color = "red")
    twin3.set_yticks(labels = 0:100:250, ticks = 0:100:250, fontsize = 12, color = "red")
end
tight_layout()

savefig(plotsdir("EGWG" * suffix *".png"))

#=
## ============ FIGURE 7: COMPARISON WITH EN4 MLD T and S ==== #


en4mld_time = 1900:2022

anom_ind = findall(x->1899<x<1970, en4mld_time)
mldT = mean(mldtemp[:, :, anom_ind], dims = 3)[:, :, 1]
mldS = mean(mldsalinity[:, :, anom_ind], dims = 3)[:, :, 1]

inds = 1900yr:10yr:1970yr


proj = ccrs.Orthographic(central_longitude=-80+55,central_latitude = 35)
noproj = ccrs.PlateCarree(central_longitude = 0)
#(80,35,-80,30,5)
lon_box = [-80, 30]
lat_box = [35, 80]

Tu =  Array(allc.ũ.dims[1]) 
Ty = Array(allc.y.dims[1])
anom_index = findall(x->1900yr<x<1970yr, Tu)

for v in [:θ, :δ]
    fig, axs = subplots(nrows = 2, ncols = 4,
                    subplot_kw = Dict("projection"=>proj), figsize = (15, 8))

    lev = v == :θ ? (-1.5:0.2:1.5) : (-1.5:0.2:1.5) ./5
θ_ = [ustrip.(value.(x)) for x in getproperty(allc, v)]
θ = [x .- (sum(θ_[anom_index]) ./ length(anom_index)) for x in θ_]
            
    for (ii, y) in enumerate(inds)
        index = findall(x->x == y, Tu)[1]
        plotme = γreshape(θ[index], allc.γ)' 
        i = ii > 4 ? 2 : 1
        j = ii > 4 ? ii - 4 : ii 
        ax =  axs[i, j]
        ax_hdl = ax.plot(orthographic_axes(lat_box..., lon_box..., 5)...,
                         color="black", linewidth=0.5,
                         transform=noproj)
        tx_path = ax_hdl[1]._get_transformed_path()
        path_in_data_coords, _ = tx_path.get_transformed_path_and_affine()
        polygon1s = mpath.Path(path_in_data_coords.vertices)
        ax.set_boundary(polygon1s)
        ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "110m", edgecolor="k", facecolor="gray"))
        
        ax.set_title(string(convert(Int64, ustrip(Tu[index]))))
        
            #add a cyclic point
        lon = ustrip.(allc.γ.lon)
        lat = ustrip.(allc.γ.lat)
        lat_box_= [lat_box[1] - 20, lat_box[2] + 20]
        lon_box_ = [360-80 - 20, 30 + 20]
        plotme = geospatialsubset(plotme', lat, lon, lat_box_, lon_box_)
        plotme = plotme'
        if v == :δ
            plotme ./= 0.5255
        end
        
        cf = ax.contourf(lon, lat, plotme, cmap = cm.balance, levels = lev, transform = noproj)
        c = ax.contour(lon, lat, plotme, colors = "black", levels = lev, transform = noproj)
        ax.clabel(c)

        #scatter
        tind = findall(x->ustrip(y)< x<ustrip(y) + 10, year.(recstime))
        recs_of_interest = vcat([r.NATL for r in recs[tind]]...)
        
        lons = Vector{Float32}(undef, 0)
        lats = Vector{Float32}(undef, 0)
        anomT = Vector{Float32}(undef, 0)
        anomS = Vector{Float32}(undef, 0)
        for r in recs_of_interest
            if !isnan(r.MLDT)
                lonind = findclosest(r.lon, lon_rs)
                latind = findclosest(r.lat, lat_rs)
                push!(lons, r.lon)
                push!(lats, r.lat)
                push!(anomT, r.MLDT + 273.15 - mldT[lonind, latind])
                push!(anomS, r.MLDS  - mldS[lonind, latind])
            end
        end
        

        stepsize = 10
        binxleft = lon_rs[begin:stepsize:end-stepsize]
        binxright = lon_rs[stepsize:stepsize:end]
        binybottom = lat_rs[begin:stepsize:end-stepsize]
        binytop = lat_rs[stepsize:stepsize:end] 
        binxcenter = (binxleft .+ binxright) ./ 2
        binycenter = (binybottom .+ binytop) ./ 2
        binned_vals = Matrix{Float32}(undef, length(binxcenter), length(binycenter))
        for (i, (lonleft, lonright)) in enumerate(zip(binxleft, binxright))
            for (j, (latbottom, lattop)) in enumerate(zip(binybottom, binytop))
                ind = intersect(findall(x->lonleft < x < lonright, lons),
                                findall(x->latbottom < x < lattop, lats))
                if v == :θ
                    binned_vals[i, j] = length(ind) > 0 ? NaNMath.mean(anomT[ind]) : NaN
                else
                    binned_vals[i, j] = length(ind) > 0 ? NaNMath.mean(anomS[ind]) : NaN
                end
                
            end
        end
        
        lons_scatter = repeat(binxcenter, outer = (1,length(binycenter)))
lats_scatter = repeat(reshape(binycenter, (1, length(binycenter))), outer =  (length(binxcenter), 1))
        nandices = findall(x->!isnan(x), binned_vals[:]) 
        s = ax.scatter(lons_scatter[nandices], lats_scatter[nandices], c = binned_vals[nandices], cmap = cm.balance, vmin = minimum(lev), vmax = maximum(lev), transform = noproj, s = 40, edgecolors = "black")
        gl = ax.gridlines(draw_labels = true)
        gl.top_labels = false
        gl.right_labels = false
        gl.bottom_labels = false
        
    end
tight_layout()
    savefig(plotsdir("post1900" * suffix * string(v) * ".png"))
end
=#



#=
d = datadir("b.e11.BLMTRC5CN.f19_g16.VOLC_GRA.001.pop.h.SST.")
cesmlme, cesmlon, cesmlat = loadCESMLME(d, 1100, 1170)
cesmlme[ismissing.(cesmlme)] .= NaN
cesmlme = convert(Matrix{Float32}, cesmlme)
figure();scatter(cesmlon, cesmlat, c = cesmlme)
=#
