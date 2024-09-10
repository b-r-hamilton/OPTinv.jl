import Pkg
Pkg.activate("../")

using OPTinv 
using Unitful, DimensionalData, BLUEs, Statistics, UnitfulLinearAlgebra, LinearAlgebra, Revise, DrWatson, TMI, CSV, DataFrames, JLD2, NaNMath, Dates, RollingFunctions, PaleoData, DateFormats, Measurements, GH19, PythonPlotExt, PythonPlot
import Measurements.value as value
import OPTinv.Est

allcores = Symbol.("MC" .* string.([28, 26, 25, 22, 21, 20, 19, 10, 9,13,14]) .* "A")
oldcores = Symbol.("MC" .* string.([26, 22, 13]) .* "A")

shallow = Symbol.("MC" .* string.([28, 26, 25]) .* "A")
intermediate = Symbol.("MC" .* string.([22, 21, 20, 19]) .* "A")
deep = Symbol.("MC" .* string.([10, 9,13,14]) .* "A")

@time allc = invert(inversion(allcores, "svd", "all", "red"))
@time oldc = invert(inversion(oldcores, "svd", "old", "blue"))

solutions = [oldc, allc]
suffix = "allvold"
#=
# ============ FIGURE 1:SOLUTIONS RECONSTRUCTED AT CORE SITES ===================== #
locs = core_locations() 
min_age = [minimum(s.y.dims[1]) for s in solutions]
max_age = [maximum(s.y.dims[1]) for s in solutions]

xl = [minimum(min_age), maximum(max_age)]

#what does the "old" solution look like, propagated to all core sites? 
ỹold = allc.predict(value.(oldc.ũ.x))
Cññold =allc.E*oldc.ũ.C*allc.E'
ỹold = DimEstimate(allc.E*oldc.ũ.v, parent(Cññold) .* permil^2, allc.y.dims)
n = UnitfulMatrix(allc.y.v) - ỹold.v
Cnn = UnitfulMatrix(ustrip.(allc.y.C), fill(permil, length(n)), fill(permil^-1, length(n)))
Jdata = n' * inv(Cnn) * n
rmserror = n' * n

fig = figure(figsize = (8, 10))
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
            #plot(s.ỹ₀.x[:, At(core)], label = "ỹ₀", color = "green", zorder = j)
        end
        if core in s.y.dims[2]
            plot(s.ỹ.x[:, At(core)], label = "ỹ: " * s.name, color = s.color, zorder = j)
        end
        if core ∉ s.y.dims[2]
            plot(ỹold.x[:, At(core)], color = "blue", zorder = j, linestyle = "dashed")
        end 
    end
    ylim(-0.2, 0.2)
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
        ylabel(L"\mathrm{\delta}^{18}\mathrm{O}_\mathrm{calcite}" * " [‰]", fontsize = 15)
    else
        ylabel("")
        yticks(-0.2:0.2:0.2, ["", "", ""])
    end
        
    title(string(core) * ", " * string(locs[core][3]) * "m", fontsize = 15)
    #grid(zorder = 1)
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

        ax = gca()
        twin2 = ax.twinx()
        twin2.spines.right.set_position(("axes", 1.15))
        steinhilber = loadSteinhilber2009()
        st_time = 1950 .- steinhilber[!, "YearBP"]
        st = Measurements.measurement.(steinhilber[!, "dTSI"], steinhilber[!, "dTSI_sigma"])
        twin2.plot(st_time, steinhilber[!, "dTSI"], color = "grey")
        twin2.fill_between(x = st_time, y1 = steinhilber[!, "dTSI"] .-  steinhilber[!, "dTSI_sigma"], y2 = steinhilber[!, "dTSI"] .+  steinhilber[!, "dTSI_sigma"], color = "grey", alpha = 0.2)
        twin2.set_ylabel("Total Solar Insolation [Wm⁻²]", fontsize = 15, color = "grey")
        twin2.set_yticks(labels = -1.5:0.5:1.5, ticks = -1.5:0.5:1.5, fontsize = 12, color = "grey")
    end
    tight_layout()
    savefig(plotsdir("modemags"*string(s)* suffix* ".png"))
end


# ================= FIGURE 3: REGION MEAN ====================== #
Tm1 = ustrip.(minimum([sol.y.dims[1][1] for sol in solutions]))
Tm2 = ustrip.(maximum([sol.y.dims[1][end-1] for sol in solutions]))
regionmeanindices = γbox(oldc.γ, 49, 89, 309, 21)
figure()
#grid()
mat = vcat(fill(1, 11113)', fill(0, (10, 11113)))
mat = vcat(fill(0,(3, 11113)), fill(1, (1, 11113)), fill(0, (7,11113)))

for sol in solutions
    sol_anom_ind = findall(x->1850yr<x<1970yr, Array(sol.ũ.dims[1]))
    θbox = estimate(sol.ũ, sol.spatialmodes, :θ, spatialinds = regionmeanindices, rolling = 2)
    θbox_alt = estimate(sol.ũ, sol.spatialmodes .* mat, :θ, spatialinds = regionmeanindices)
    μ = mean(θbox.v[sol_anom_ind])
    θbox = DimEstimate(θbox.v .- μ, θbox.C, θbox.dims)
    θbox_alt = DimEstimate(θbox_alt.v .- mean(θbox_alt.v[sol_anom_ind]), θbox_alt.C, θbox_alt.dims)

    inds = findall(x->x∈sol.y.dims[1], Array(sol.ũ.dims[1]))
    pre1900inds = intersect(inds, findall(x->1000yr<x<1900yr, Array(sol.ũ.dims[1])))
    minv = findmin(θbox.x[inds])
    println(sol.name * ": minima")
    @show sol.ũ.dims[1][inds][minv[2]]
    
    println(sol.name * ": pre-1900 maxima") 
    @show max = findmax(θbox.x[pre1900inds])
    @show sol.ũ.dims[1][pre1900inds][max[2]]

    println("LLS, LIA cooling")
    #we need to subset to the time period we want to compute the slope over 
    liainds = findall(x->sol.ũ.dims[1][inds][minv[2]] > x > sol.ũ.dims[1][pre1900inds][max[2]], Array(sol.ũ.dims[1]))
    x = UnitfulMatrix(Array(sol.ũ.dims[1][liainds]))
    C = UnitfulMatrix(parent(θbox.C[liainds, liainds]), unitrange(θbox.C)[liainds], unitdomain(θbox.C)[liainds])
    y = UnitfulMatrix(θbox.v[liainds])
    lls, C = linearleastsquares(x, y, C=C)
    @show lls[1] * 1000 
    @show sqrt.(diag(C))[1] .* 1000

    plot(θbox.x[inds], color = sol.color, label = sol.name, zorder = 10000)
    newcolor = sol.color == "red" ? "pink" : "aqua"
    #plot(value.(θbox_alt.x[inds]), color = newcolor, label = sol.name * ", mode 1 recon.", linewidth = 3, zorder = 10000, linestyle = "dashed")
    println()
end

t_ocean2k = 100:200:1900
oc2k_binned, oc2k_binnedv, oc2k_ages = loadOcean2kBinned()
#the following values are from from McGregor 2015, Supp. Table S13
#sd_ocean2k = vec([0.78 0.58 0.39 0.38 0.23 0.07 -0.19 -0.70 -0.71 -0.60])
#this is the same as the above (and I promise the timing matches up) 
median_sd_ocean2k = reverse([NaNMath.median(oc2k_binned[i, :]) for i in 1:size(oc2k_binned)[1]])
#Here, I calculate the 0.25 and 0.75 percentile range (same as McGregor 2015 Figure 2)
std_sd_ocean2k = reverse([quantile(oc2k_binned[i,:][(!).(isnan.(oc2k_binned[i, :]))], [0.25, 0.75]) for i in 1:size(oc2k_binned)[1]])

#the following 0.42degC/kyr value is from from McGregor 2015, Supp. Table S6
#this calculation converts to °C 
#0.42degC/kyr * (sd) / (std/kyr) 
y_ocean2k = 0.42.*median_sd_ocean2k./(median_sd_ocean2k[4]-median_sd_ocean2k[9])
#apply same conversion to the quantile values 
ystd_ocean2k = 0.42 .* std_sd_ocean2k ./(median_sd_ocean2k[4]-median_sd_ocean2k[9])
ystd_ocean2k = hcat(ystd_ocean2k...) #2 x 10
#to plot these, we need them as +/- from the median estimate values
ystd_ocean2k[1, :] .= y_ocean2k .- ystd_ocean2k[1, :]
ystd_ocean2k[2, :] .=  ystd_ocean2k[2, :] .- y_ocean2k
#we want this to be anomaly from 1850-1970, but we can't do that for this dataset, next best thing is just use the last bin 
y_ocean2k .-= y_ocean2k[end]
#then use some hlines bb 
hlines(y = y_ocean2k, xmin = t_ocean2k .- 100, xmax = t_ocean2k .+ 100, color = "gray", linewidth = 3)
hlines(y = y_ocean2k .+ ystd_ocean2k[2, :], xmin = t_ocean2k .- 100, xmax = t_ocean2k .+ 100, color = "gray")
hlines(y = y_ocean2k .- ystd_ocean2k[1, :], xmin = t_ocean2k .- 100, xmax = t_ocean2k .+ 100, color = "gray")
vlines(x = t_ocean2k .- 100, ymin = y_ocean2k .- ystd_ocean2k[1, :], ymax = y_ocean2k .+ ystd_ocean2k[2,:], color = "gray")
vlines(x = t_ocean2k .+ 100, ymin = y_ocean2k .- ystd_ocean2k[1, :], ymax = y_ocean2k .+ ystd_ocean2k[2,:], color = "gray")
[fill_between(x = [t-100, t+100], y1 = y .- ystd1, y2 = y .+ ystd2, color = "gray", alpha = 0.25) for (t, y, ystd1, ystd2) in zip(t_ocean2k, y_ocean2k, ystd_ocean2k[1, :], ystd_ocean2k[2, :])]

#errorbar(x = t_ocean2k, y = y_ocean2k, yerr = ystd_ocean2k)

#GH19
exp = explist()[end]
filename = download_exp(exp, true)
using NCDatasets
tmi = NCDataset(TMI.download_ncfile(TMIversion()))
lat = tmi["lat"][:]
lon = tmi["lon"][:]
natl = tmi["d_NATL"][:, :]
natl[natl .== 0] .= NaN
ant =  tmi["d_ANT"][:, :]
ant[ant .== 0] .= NaN
gin = tmi["d_GIN"][:,:]
gin[gin .== 0] .= NaN               
nc = NCDataset(filename)
θGH19 = nc["theta"][:, 1, :, :]
tGH19 = nc["year"][:]
θGH19NATL = [NaNMath.mean(θGH19[i,:, :] .* natl') for i in 1:size(θGH19)[1]]
θGH19ANT = [NaNMath.mean(θGH19[i,:, :] .* ant') for i in 1:size(θGH19)[1]]
θGH19GIN = [NaNMath.mean(θGH19[i,:, :] .* gin') for i in 1:size(θGH19)[1]]
ind1970 = findall(x->x==1970, tGH19)[1]
ind1880 = findall(x->x==1880, tGH19)[1]
ind1780 = findall(x->x==1780, tGH19)[1]    
θGH19NATL[ind1970]-θGH19NATL[ind1780]
θGH19NATL[ind1970]-θGH19NATL[ind1780]

#figure();plot(tGH19, θGH19NATL, color = "pink", label = "NATL");plot(tGH19, θGH19ANT, color = "aqua", label = "ANT");plot(tGH19, θGH19GIN, color = "green", label = "GIN");legend()
#title("Mean Temperature Anomaly from 1850-1970yr: 50°W-20°E, 50°N-90°N") 
lmrdataset = loadLMR("sst")
lmrtime = lmrdataset["time"][:]; lmrlon = lmrdataset["lon"][:]; lmrlat = lmrdataset["lat"][:]

for i in 1:20
    lmrsst = makeNaN(lmrdataset["sst"][:, :, i, :])
    anomindex = findall(x->1850 < x < 1970, year.(lmrtime))
    lmrsst = lmrsst .- mean(lmrsst[:, :, anomindex], dims = 3)[:, :, 1]

    lmrtimeroll = rolling(mean, year.(lmrtime), 50)
    lmrglobal = rolling(mean, [NaNMath.mean(lmrsst[:, :, i]) for i in 1:size(lmrsst)[3]], 50)

    lmrbox = rolling(mean, boxmean(lmrsst, lmrlat, lmrlon, 50, 90, 360-50, 20), 50)
    #plot(DimArray((lmrglobal .- mean(lmrglobal)) * K, Ti(lmrtime)), color = "black", label = "LMR: global", zorder = 0)
    plot(DimArray((lmrbox  .- mean(lmrbox))* K, Ti(lmrtimeroll)), color = "black", label = "LMR", zorder = 0, linewidth = 1, linestyle = "solid")
end
lmrsst = makeNaN(mean(lmrdataset["sst"][:, :, :, :], dims = 3))[:, :, 1, :]
    
hadisst = loadHadISST()
hadsst = makeNaN(hadisst["tos"][:, :, :])
hadlat = hadisst["latitude"][:]; hadlon = hadisst["longitude"][:]; hadtime = hadisst["time"][:]
hadroll(x) = rolling(mean, x, 50 * 12)
hadisst_gm = hadroll([NaNMath.mean(hadsst[:, :, i]) for i in 1:size(hadsst)[3]])
hadisst_bm = hadroll(boxmean(hadsst, hadlat, hadlon, 50, 90, -50, 20))
hadtimeroll = hadroll(yeardecimal.(hadtime))
#plot(DimArray(hadisst_gm .- hadisst_gm[begin], Ti(hadtimeroll)), color = "green", label = "HadISST: global")
#plot(DimArray(hadisst_bm .- mean(hadisst_bm), Ti(hadtimeroll)), color = "black", label = "HadISST", linewidth = 3, linestyle = "dashed")

ylabel("Surface Temperature Anomaly [K]", fontsize = 15)
xlabel("Time [yr CE]", fontsize = 15)
#=
oc2k_binned, oc2k_binnedv, oc2k_ages = loadOcean2kBinned()
oc2k, oc2kdata = loadOcean2k()

atl_indices = findall(x->occursin("Atlantic", x), oc2k[!, "name"])
arc_indices = findall(x->occursin("Arctic", x), oc2k[!, "name"])
inregion = findall(x->x>50, oc2k[!, "lat"])
rel_records = intersect(vcat(atl_indices, arc_indices), inregion)

oc2kglobal = [NaNMath.mean(oc2k_binned[i, :]) for i in 1:size(oc2k_binned)[1]]
oc2kglobalstd = [NaNMath.std(oc2k_binned[i, :]) for i in 1:size(oc2k_binned)[1]]
oc2kregion = [NaNMath.mean(oc2k_binned[i, rel_records]) for i in 1:size(oc2k_binned)[1]]




figure();scatter(reverse(collect(bincenter)), Measurements.measurement.(oc2kglobal, oc2kglobalstd));ylim(-3,3)
oc2kglobal
lls, Cnew = linearleastsquares(reverse(collect(bincenter)), oc2kglobal, C = Diagonal(oc2kglobalstd)) #nice!
@show Measurements.measurement(lls[2],sqrt(Cnew[2,2]))



binleft = 0:200:2000
bincenter = 100:200:2000
binright = 200:200:2000
size(oc2k)
oc2kanom = Matrix{Any}(undef, size(oc2k)[1], length(bincenter))
oc2kanom .= NaN
#oc2kanomstd = copy(oc2kanom)

for i in atl_indices
    sstname = names(oc2kdata)[i]
    meta_index = findall(x->occursin(sstname[begin:begin+11], x), oc2k.name)[1]
     
    t = oc2kdata[!, sstname]
    t[t .== "NaN"] .= NaN
    y = oc2kdata[!, names(oc2kdata)[i*2]]
    y[y .== "NaN"] .= NaN
    for (j, (bl, br)) in enumerate(zip(binleft, binright))
        ind = findall(x->bl<x<br, t)
        if length(ind) > 0
            oc2kanom[i, j] = Measurements.measurement(mean(y[ind]), std(y[ind]))
            
        else
            oc2kanom[i, j] = NaN
        end
    end
end
oc2kanom[(!).(isnan.(oc2kanom[:, 1]))]
[NaNMath.mean(oc2kanom[:, i]) for i in 1:10]

#oc2kanom .-= oc2kanom[:, end]

#figure()

#[plot(bincenter, oc2kanom[i, :], color = "gray") for i in 1:57]
#plot(bincenter, [NaNMath.mean(oc2kanom[:, i]) for i in 1:10], color = "black", label = "Ocean2k: Atlantic Mean", ".-")

=#

#=
slope_est = (0.47, 0.36)
t = Array(ustrip.(oldc.ũ.dims[1]))
sol_anom_ind = findall(x->1850<x<1970, t)
t_oc2k = 1500yr:1yr:1700yr
y1 = [(1970 - t_) * slope_est[1] * 0.001 for t_ in ustrip.(t_oc2k)]
y2 = [(1970 - t_) * slope_est[2] * 0.001 for t_ in ustrip.(t_oc2k)]
#fill_between(x = t, y1 = y1 .- mean(y1[sol_anom_ind]), y2 = y2 .- mean(y2[sol_anom_ind]), color = "black", alpha = 0.5, label = "Global Ocean2k SST Trend")

plot(ustrip.(t_oc2k), y1 .+ 0.05, color = "black", label = "Global Oc2k SST Trend")
text(x = 1709, y = 0.135, s = "0.47K/kyr:global")
plot(ustrip.(t_oc2k), y2 .+ 0.05, color = "black")
text(x = 1709, y = 0.18, s = "0.36K/kyr:N.Atl.")
=#
#legend()
xlim(Tm1, Tm2)
xticks(800:200:1800, fontsize = 12) 
ylim(-0.5, 0.9)
yticks(-0.5:0.25:0.75, fontsize = 12)
tight_layout()
savefig(plotsdir("meants" * suffix * ".png"))

# ================= FIGURE 4: MAPS ============= #
inds = collect(1100yr:200yr:2000yr)
stepsize = 70  
lev = -1.5:0.2:1.5

proj = ccrs.Orthographic(central_longitude=-80+55,central_latitude = 35)
noproj = ccrs.PlateCarree(central_longitude = 0)
fig, axs = subplots(nrows = length(solutions), ncols = length(inds),
                    subplot_kw = Dict("projection"=>proj), figsize = (15, 4),
                    constrained_layout = true)
#(80,35,-80,30,5)
lon_box = [-80, 30]
lat_box = [35, 80]
#for (i, sol) in enumerate([solutions..., lmrsst])
for (i, sol) in enumerate(solutions) 
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
        #PyPlot version
        #ax =  axs[(ii-1) * (length(solutions)+1)+ i]
        #ax =  axs[(ii-1) * (length(solutions))+ i]
        #PythonCall version
        ax = axs[i-1, ii-1]
        if (sol isa solution && y ∈ Ty) || (sol isa Array)
            
            ax_hdl = ax.plot(orthographic_axes(lat_box..., lon_box..., 5)...,
                             color="black", linewidth=0.5,
                             transform=noproj)
            #PyPlot version
            #tx_path = ax_hdl[1]._get_transformed_path()
            #PythonPlot version
            tx_path = first(ax_hdl)._get_transformed_path()
            
            path_in_data_coords, _ = tx_path.get_transformed_path_and_affine()
            polygon1s = mpath.Path(path_in_data_coords.vertices)
            ax.set_boundary(polygon1s)
            ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "110m", edgecolor="k", facecolor="gray"))

            if i == 1
                ax.set_title(string(convert(Int64, ustrip(Tu[start_index]))) * "-" *string(convert(Int64, ustrip(Tu[stop_index]))), fontsize = 20, fontweight = "bold")
            end
            ax.set_extent([lon_box..., lat_box[1]- 10, lat_box[2]], noproj)
            
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
        else
            ax.clear()
            ax.axis("off") 
        end
    end
    tight_layout()
    savefig(plotsdir("surfacesol" * suffix * ".png"), bbox_inches = "tight")
end

#=
# ==================== Figure 4alt: MAP of MCA-LIA  ==================== #
lev = -1.5:0.2:1.5

proj = ccrs.Orthographic(central_longitude=-80+55,central_latitude = 35)
noproj = ccrs.PlateCarree(central_longitude = 0)
fig, axs = subplots(nrows = 1, ncols = 1,
                    subplot_kw = Dict("projection"=>proj), figsize = (10, 6))
#(80,35,-80,30,5)
lon_box = [-80, 30]
lat_box = [35, 80]

θ_ = [ustrip.(value.(x)) for x in oldc.θ]
mcaind = findall(x->950yr < x < 1250yr, Array(oldc.ũ.dims[1]))
liaind = findall(x->1450yr < x < 1850yr, Array(oldc.ũ.dims[1]))

plotme = γreshape(mean(θ_[mcaind]) - mean(θ_[liaind]), oldc.γ)' #100yr at 10yr res

ax =  axs
ax_hdl = ax.plot(orthographic_axes(lat_box..., lon_box..., 5)...,
                 color="black", linewidth=0.5,
                 transform=noproj)
tx_path = ax_hdl[1]._get_transformed_path()
path_in_data_coords, _ = tx_path.get_transformed_path_and_affine()
polygon1s = mpath.Path(path_in_data_coords.vertices)
ax.set_boundary(polygon1s)
ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "110m", edgecolor="k", facecolor="gray"))
            
#ax.set_title(string(convert(Int64, ustrip(Tu[start_index]))) * "-" *string(convert(Int64, ustrip(Tu[stop_index]))))
            
#add a cyclic point
lat_box_= [lat_box[1] - 20, lat_box[2] + 20]
lon_box_ = [360-80 - 20, 30 + 20]

plotme = geospatialsubset(plotme', oldc.γ.lat, oldc.γ.lon, lat_box_, lon_box_)
plotme = plotme'
cf = ax.contourf(oldc.γ.lon, oldc.γ.lat, plotme, cmap = cm.balance, levels = lev, transform = noproj)
c = ax.contour(oldc.γ.lon, oldc.γ.lat, plotme, colors = "black", levels = lev, transform = noproj)
ax.clabel(c)
gl = ax.gridlines(draw_labels = true)
gl.top_labels = false
gl.right_labels = false
gl.bottom_labels = false


tight_layout()
savefig(plotsdir("mcalia.png"))

# ==================== Figure 4alt2: MAP of PD-LIA  ==================== #
lev = -1.5:0.2:1.5

proj = ccrs.Orthographic(central_longitude=-80+55,central_latitude = 35)
noproj = ccrs.PlateCarree(central_longitude = 0)
fig, axs = subplots(nrows = 1, ncols = 2,
                    subplot_kw = Dict("projection"=>proj), figsize = (15, 8))
#(80,35,-80,30,5)
lon_box = [-80, 30]
lat_box = [35, 80]
for (i, sol) in enumerate([oldc, allc])
    θ_ = [ustrip.(value.(x)) for x in sol.θ]
    mcaind = findall(x->1900yr < x < 1950yr, Array(sol.ũ.dims[1]))
    liaind = findall(x->1500yr < x < 1850yr, Array(sol.ũ.dims[1]))

    plotme = γreshape(mean(θ_[mcaind]) - mean(θ_[liaind]), sol.γ)' #100yr at 10yr res

    ax =  axs[i]
    ax_hdl = ax.plot(orthographic_axes(lat_box..., lon_box..., 5)...,
                     color="black", linewidth=0.5,
                     transform=noproj)
    tx_path = ax_hdl[1]._get_transformed_path()
    path_in_data_coords, _ = tx_path.get_transformed_path_and_affine()
    polygon1s = mpath.Path(path_in_data_coords.vertices)
    ax.set_boundary(polygon1s)
    ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "110m", edgecolor="k", facecolor="gray"))
    
    #ax.set_title(string(convert(Int64, ustrip(Tu[start_index]))) * "-" *string(convert(Int64, ustrip(Tu[stop_index]))))
    
    #add a cyclic point
    lat_box_= [lat_box[1] - 20, lat_box[2] + 20]
    lon_box_ = [360-80 - 20, 30 + 20]

plotme = geospatialsubset(plotme', sol.γ.lat, sol.γ.lon, lat_box_, lon_box_)
plotme = plotme'
cf = ax.contourf(sol.γ.lon, sol.γ.lat, plotme, cmap = cm.balance, levels = lev, transform = noproj)
c = ax.contour(sol.γ.lon, sol.γ.lat, plotme, colors = "black", levels = lev, transform = noproj)
ax.clabel(c)
gl = ax.gridlines(draw_labels = true)
gl.top_labels = false
gl.right_labels = false
gl.bottom_labels = false


    tight_layout()
end

savefig(plotsdir("pdlia.png"))
=#
# ===================== FIGURE 5alt: Labrador T and S  ================ #

df = CSV.read(OPTinv.datadir("TLS1998Fig7Digitized.csv"), DataFrame)

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

df = CSV.read(OPTinv.datadir("TLS1998Fig7Digitized.csv"), DataFrame)
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
rolling(NaNMath.mean, hadsst_lab, 10)

plot( hadtimeann, hadsst_lab, color = "black", label = "HadiSST", zorder = 1, linestyle = "dashed")
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

fig = figure(figsize = (10,8));
fig.subplots_adjust(right=0.75)
for (i, sol) in enumerate(solutions)
    
    ax = subplot(length(solutions), 1, i)
    grid()
    T = sol.y.dims[1]
    println("SOLUTION = " * sol.name)
    for (k, ls, lw, col) in zip(keys(inds), ["solid", "solid", "dotted", "dashed"], [4,4,3,3], ["purple", "darkgreen", "darkgreen", "darkgreen"])
        Tu = sol.ũ.dims[1]
        temp = estimate(sol.ũ, sol.spatialmodes, :θ, spatialinds = inds[k])
        anom = mean(ustrip.(temp.x[DimensionalData.Between(1850yr, 1970yr)]))
        y1 = ustrip.(temp.x[DimensionalData.Between(T[1], T[end])]) .- anom 
        plot(ustrip.(T), value.(y1), color = col, label = k, linestyle = ls, linewidth = lw)
        
        println(string(k))
        @show min = findmin(y1)
        @show t1 = T[min[2]]
        pre1900inds = findall(x->x<1900yr, Array(T))
        @show max = findmax(y1[pre1900inds])
        @show t2 = T[pre1900inds][max[2]]
        println(estimate(sol.ũ, sol.spatialmodes, :θ, t1,t2,Tm1, Tm2, spatialinds = inds[k]))
        @show t1
        t1ind = findall(x->x == t1, Array(Tu))[1]
        t2ind = findall(x->x == t2, Array(Tu))[1]
        #calculate slope 
        wind = t2ind:t1ind
        x = UnitfulMatrix(Array(sol.ũ.dims[1][wind]))
        C = UnitfulMatrix(parent(temp.C[wind, wind]), unitrange(temp.C)[wind], unitdomain(temp.C)[wind])
        y = UnitfulMatrix(temp.v[wind])
        
        lls, C = linearleastsquares(x, y; C)
        @show getindexqty(lls,1) * 1000 
        @show sqrt.(diag(C))[1] .* 1000

        println()
        
    end
    # plot(tGH19, θGH19NATL, color = "pink")
    # plot(tGH19, θGH19ANT, color = "aqua")
    #title(sol.name)
    xlim(Tm1, Tm2)
    #legend(loc = "upper right")
    ylabel("Temperature Anomaly [K]",fontsize = 15)
    if i == 2 xlabel("Time [years CE]", fontsize = 15) end
 
    yticks(-0.5:0.25:0.75, fontsize = 12)
    ylim(-0.5, 0.8)
    xticks(800:200:1800, fontsize = 12)
    text(x = 810, y = -0.48, s = ["A", "B"][i], fontsize = 30, weight = "bold")
    #=
    twin1 = gca().twinx()
    dahljensen = CSV.read(OPTinv.datadir("DigitizedDahlJensen.csv"), DataFrame)
    djyears = dahljensen[!, "Time [years CE]"]
    anom_index = findall(x->1850 < x < 1970, djyears)
    djtemp = dahljensen[!, "Temperature [degC]"] #.- mean(dahljensen[anom_index, "Temperature [degC]"])  
    twin1.plot(djyears, djtemp, color = "blue", linewidth = 1, alpha = 0.5)
    twin1.set_ylabel("Borehole Temp. [°C]", fontsize = 15, color = "blue")
    twin1.set_yticks(labels = -32.5:0.5:-30, ticks = -32.5:0.5:-30, fontsize = 12, color = "blue")
    =#
    twin2 = ax.twinx()
    #twin2.spines.right.set_position(("axes", 1.15))
    steinhilber = loadSteinhilber2009()
    st_time = 1950 .- steinhilber[!, "YearBP"]
    st = Measurements.measurement.(steinhilber[!, "dTSI"], steinhilber[!, "dTSI_sigma"])
    twin2.plot(st_time, steinhilber[!, "dTSI"], color = "grey")
    twin2.fill_between(x = st_time, y1 = steinhilber[!, "dTSI"] .-  steinhilber[!, "dTSI_sigma"], y2 = steinhilber[!, "dTSI"] .+  steinhilber[!, "dTSI_sigma"], color = "grey", alpha = 0.2)
    twin2.set_ylabel("Total Solar Insolation [Wm⁻²]", fontsize = 15, color = "grey")
    twin2.set_yticks(labels = -1.5:0.5:1.5, ticks = -1.5:0.5:1.5, fontsize = 12, color = "grey")

    twin3 = ax.twinx()
    twin3.spines.right.set_position(("axes", 1.1))
    gao = loadGao2008()
    twin3.plot(gao[!, "Year"], gao[!, "Global"], color = "red", alpha = 0.5, linewidth = 1)
    twin3.set_ylabel("Global strat. sulf. aer. inj. [Tg]", fontsize = 15, color = "red")
    twin3.set_yticks(labels = 0:100:250, ticks = 0:100:250, fontsize = 12, color = "red")
end
tight_layout()
savefig(plotsdir("EGWG" * suffix *".png"))


## ============ X : EFFECTIVE δ¹⁸Oc ============= ##
modes = allc.spatialmodes
cdims = vec(covariancedims(allc.ũ.dims))
Tu = Array(allc.ũ.dims[1])
Vsub = UnitfulMatrix(modes[:, regionmeanindices], fill(NoUnits, 11), fill(NoUnits, 500))
Ey = fill(1/500, 500)
Edag = Vsub * Ey

mat = Matrix{Float64}(undef,length(Tu),length(cdims))
mat .= 0
for (i,t) in enumerate(Tu)
    xindsθ = findall(x->x[1]==t && x[3] == :θ,cdims)
    xindsδ = findall(x->x[1]==t && x[3] == :δ,cdims)
    mat[i, xindsθ] = Edag * -0.224
    mat[i, xindsδ] = Edag
end

mat = UnitfulMatrix(mat, fill(permil, length(Tu)), unitrange(allc.ũ.v))
figure()
for (tit, sol) in zip(["old", "all"], solutions)
    Css = mat*sol.ũ.C*transpose(mat)
    y = UnitfulMatrix(vec(mat*sol.ũ.v))
    unc = UnitfulMatrix(sqrt.(diag(parent(Css))), unitrange(y))
    plot(DimArray(measurement.(vec(y), vec(unc)), Ti(Array(sol.ũ.dims[1]))), label = tit, color = sol.color)
    #legend()
    if tit == "old"
        
        ind = findall(x->1900yr>x>1140yr, Tu)
        display(parent(y)[ind])
        l, c= linearleastsquares(ustrip.(Array(Tu))[ind], parent(y)[ind], C = parent(Css)[ind,ind])
        @show l
        display(sqrt.(diag(c)))
        
    end
end

xlabel("Time [years CE]", fontsize = 15)
ylabel("Effective " * L"\mathrm{\delta}^{18}\mathrm{O_{calcite}}" * " [‰]", fontsize = 15)
xticks(600:200:2000, fontsize = 12)
yticks(-0.2:0.1:0.2, fontsize = 12)
gca().invert_yaxis()
tight_layout()
savefig(plotsdir("effd18Omean.png"))

# ============ FIGURE X: PLANKTIC STACK ================= #
cdims = vec(covariancedims(allc.ũ.dims))
Tu = Array(allc.ũ.dims[1])
corelons = 360 .+ [c[1] for c in locs]
corelats = [c[2] for c in locs]
surfind = γbox(oldc.γ, minimum(corelats)-1, maximum(corelats)+1, minimum(corelons)-1, maximum(corelons)+1)

Vsub = UnitfulMatrix(modes[:, surfind], fill(NoUnits, 11), fill(NoUnits, 2))
Ey = fill(1/2, 2)
Edag = Vsub * Ey

mat = Matrix{Float64}(undef,length(Tu),length(cdims))
mat .= 0
for (i,t) in enumerate(Tu)
    xindsθ = findall(x->x[1]==t && x[3] == :θ,cdims)
    xindsδ = findall(x->x[1]==t && x[3] == :δ,cdims)
    mat[i, xindsθ] = Edag * -0.224
    mat[i, xindsδ] = Edag
end

mat = UnitfulMatrix(mat, fill(permil, length(Tu)), unitrange(allc.ũ.v))
figure()
for (tit, sol) in zip(["old", "all"], solutions)
    Css = mat*sol.ũ.C*transpose(mat)
    y = UnitfulMatrix(vec(mat*sol.ũ.v))
    unc = UnitfulMatrix(sqrt.(diag(parent(Css))), unitrange(y))
    plot(DimArray(measurement.(vec(y), vec(unc)), Ti(Array(sol.ũ.dims[1]))), label = tit, color = sol.color)
    if tit == "old"
        ind = findall(x->1800yr>x>1200yr, Tu)
        l, c= linearleastsquares(ustrip.(Tu)[ind], parent(y)[ind], C = parent(Css[ind,ind]))
        @show l
        display(sqrt.(diag(c)))
    end
end


#make a planktic stack
Gb = loadcores(oldcores, dir = OPTinv.datadir("runBacon"), rules = ["G.b"])
Gi = loadcores(oldcores, dir = OPTinv.datadir("runBacon"), rules = ["G.i"])

Gbunc = reshape(sqrt.(diag(parent(Gb.Cnn.mat)))permil, size(Gb.y))
Giunc = reshape(sqrt.(diag(parent(Gi.Cnn.mat)))permil, size(Gi.y))

mat = hcat(Matrix(Gi.y) .± Giunc, Matrix(Gb.y) .± Gbunc)
mat = Matrix(Gb.y) .± Gbunc
inds = findall(x->1980yr > x > 1850yr, Array(Gb.y.dims[1]))
#remove 1850-1980 mean


mat .-= mean(value.(mat[inds, :]), dims = 1)
plot(DimArray(mean(mat, dims = 2)[:], Ti(Array(Gb.y.dims[1]))), color = "black")
xlabel("Time [years CE]", fontsize = 15)
ylabel("Effective " * L"\mathrm{\delta}^{18}\mathrm{O_{calcite}}" * " [‰]", fontsize = 15)
xticks(600:200:2000, fontsize = 12)
yticks(-0.4:0.1:0.4, fontsize = 12)
gca().invert_yaxis()
tight_layout()
savefig(plotsdir("effd18Ostack.png"))

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


# ============= FIGURE:OLDC Slope at every point ====================== #
ind = findall(x->1900yr>x>1140yr, Array(Tu))
sol = oldc
θmat = value.(ustrip.(hcat(sol.θ...)))
lls = Vector{Float64}(undef, size(θmat)[1])

for i in 1:size(lls)[1]
    lls[i] = linearleastsquares(ustrip.(Array(Tu)[ind]), θmat[i, ind])[1][1]
end
proj = ccrs.Orthographic(central_longitude=-80+55,central_latitude = 35)
noproj = ccrs.PlateCarree(central_longitude = 0)
fig, axs = subplots(nrows = 1, ncols = 1,
                    subplot_kw = Dict("projection"=>proj), figsize = (10,8))
ax = axs
ax_hdl = ax.plot(orthographic_axes(lat_box..., lon_box..., 5)...,
                 color="black", linewidth=0.5,
                 transform=noproj)
tx_path = ax_hdl[0]._get_transformed_path()
path_in_data_coords, _ = tx_path.get_transformed_path_and_affine()
polygon1s = mpath.Path(path_in_data_coords.vertices)
ax.set_boundary(polygon1s)
ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "110m", edgecolor="k", facecolor="gray"))
lon = ustrip.(oldc.γ.lon)
lat = ustrip.(oldc.γ.lat)
lat_box_= [lat_box[1] - 20, lat_box[2] + 20]
lon_box_ = [360-80 - 20, 30 + 20]
plotme = γreshape(lls .* 1000, oldc.γ)
@show NaNMath.maximum(plotme)
#plotme = geospatialsubset(plotme, lat, lon, lat_box_, lon_box_)
plotme, lon_ = cu.add_cyclic_point(plotme', coord = lon)
cf = ax.contourf(lon_, lat, plotme, cmap = cm.balance, vmin = -3, vmax = 3, transform = noproj)
c = ax.contour(lon_, lat, plotme, colors = "black", vmin = -3, vmax = 3, transform = noproj) 
ax.clabel(c)

#find LLS of Ocean2K records in region of interest, if they have data between 1000 and 1500
lons = [-50, 20]
lats = [50, 90] 
oc2k_binned, oc2k_binnedv, oc2k_ages = loadOcean2kBinned()
oc2k, oc2kdata = loadOcean2k()

lon = oc2k[!, "lon"]
lat = oc2k[!, "lat"]
recinds = intersect(findall(x->-50 < x < 20, lon), findall(x->50<x<90, lat))
recs = oc2k[recinds, "name"]

function getrid(x)
    x[x.== "NaN"] .= NaN
    return convert(Vector{Float64}, x)
end

d = Dict() 
for (i, r) in enumerate(recs)
        
    locmatind = recinds[i]
    lat_ = lat[locmatind]
    lon_ = lon[locmatind]
    @show r
    inds = findall(x->occursin(r[begin:20], x), names(oc2kdata))
    cutoff = findall(x->x=='_', names(oc2kdata)[inds[1]])[1]
    lab = names(oc2kdata)[inds[1]][begin:cutoff-1]
    age = getrid(oc2kdata[!, inds[1]])
    sst = getrid(oc2kdata[!, inds[2]])
    #=
    age = collect(100.0:200:2000)
    oc2k_binned
    sst = oc2k_binned[:, locmatind]
    
    nanindex = sort(findall(x->!isnan(x), sst))
    age = age[nanindex]
    sst = sst[nanindex]
    =#
    post1000 = findall(x->1140<x<1900, age)
    if NaNMath.maximum(age) > 1500
        lls, C = linearleastsquares(age[post1000], sst[post1000])
        @show lls[1] * 1000
        d[r] = (lat_, lon_, lls[1] * 1000)
    end
end
for (i, k) in enumerate(keys(d))
    if d[k][1] < 70 #don't include disagreeing Arctic recs
        s = scatter(d[k][2], d[k][1], c = d[k][3], vmin = -3, vmax = 3, cmap = cm.balance, s = 100, edgecolors = "white", transform = noproj, zorder = 100000,linewidths = 2)
    end
end


tight_layout()
savefig(plotsdir("ocean2kcomp.png"))

# =============== MODES =============== #
figure(figsize = (10,10))
for i in 1:11
    ax = subplot(4,3,i)
    plot(allc.u₀.x[:, At(i), At(:θ)], color = "gray")
    plot(allc.ũ.x[:, At(i), At(:θ)], color = "red")
    plot(oldc.ũ.x[:, At(i), At(:θ)], color = "blue")
    title(string(i), fontsize = 15)
    yt = ax.get_yticks()
    ax.set_yticks(yt,yt,fontsize = 12)
    if i ∈ [1,4,7,10]
        ylabel("Mode Mag. [K]", fontsize = 15)
    end
    xt = collect(500:500:2000)    
    if i ∉ [9,10,11]
        xticks(xt,fill("", length(xt)))
        xlabel("")
    else
        xlabel("Time [years CE]", fontsize = 15)

        xticks(xt, xt, fontsize = 12)
    end    
end
#tight_layout()
savefig(plotsdir("u0utilde.png"))

# === just Mode 1 plot  === #
gkw = ("height_ratios" => [1,1],)
fig, axs = subplots(2,1, gridspec_kw = gkw)
function specialplot(da::DimArray, ax, color::String; fb = true)
    x = ustrip.(Array(dims(da)[1]))
    y = Array(da)
    yunc = ustrip.(Measurements.uncertainty.(y))
    y = ustrip.(value.(y))
    ax.plot(x,y, color = color)
    if fb
        ax.fill_between(x, y1 = y .- yunc, y2 = y.+yunc, color = color, alpha = 0.5)
    end
    
end

specialplot(allc.u₀.x[:, At(1), At(:θ)], axs[0], "gray")
specialplot(allc.ũ.x[:, At(1), At(:θ)], axs[0], "red")
specialplot(oldc.ũ.x[:, At(1), At(:θ)], axs[0], "blue")
axs[0].set_ylabel("Mode Mag. [K]", fontsize = 15)
yt = axs[0].get_yticks()
axs[0].set_yticks(yt,yt,fontsize = 12)
xt = collect(500:500:2000)
axs[0].set_xticks(xt, fill("", length(xt)))
axs[0].set_xlabel("")
specialplot(allc.u₀.x[:, At(1), At(:θ)], axs[1], "gray", fb = false)
specialplot(allc.ũ.x[:, At(1), At(:θ)], axs[1], "red", fb = false)
specialplot(oldc.ũ.x[:, At(1), At(:θ)], axs[1], "blue", fb = false)
axs[1].set_ylabel("Mode Mag. [K]", fontsize = 15)
axs[0].text(x = 500, y = -60, s = "A", fontsize = 30, weight = "bold")
axs[1].text(x = 500, y = -20, s = "B", fontsize = 30, weight = "bold")
xticks(xt, xt, fontsize = 12)
axs[1].set_xlabel("Time [years CE]", fontsize = 15)
tight_layout()
savefig(plotsdir("u0utilde_mode1.png"))


#=
# =============== building up SVD maps ===== #
figure(figsize = (10,10))
for i in 1:11#iterate through mode
    for j in 1:11 #iterate through cores 
        subplot(11,11, (i-1)*11+j)
        v = s.U[j,i] * s.S[j] .* s.Vt[i, :]
        pc = pcolormesh(allc.γ.lon, allc.γ.lat, γreshape(v, allc.γ)')
    end
end
tight_layout()
=#


# =============== SINGULAR VALUES =========== #
s = jldopen("../data/M/svd.jld2")["SVD"]
U = s.U #cores × modes

depths = [locs[c][3] for c in allcores]
figure(figsize = (10, 4))
for i in 1:11
    ax = subplot(1, 11, i) 
    #plot(U[i, :], Array(allc.y.dims[2]), color = cm((11-i)/11), label = i)
    scatter(U[:, i], depths, c = U[:, i], cmap = cm.balance, vmin = -1, vmax = 1,edgecolors = "gray", 100)
    ax.invert_yaxis()
    xlim(-1, 1)
    xl = ax.get_xlim()
    yt = collect(1000:200:2200)
    hlines(y = depths, xmin = xl[0], xmax = xl[1], color = "gray", zorder = 0)
    tx = ax.twinx()
    tx.set_xlim(xl)
    
    if i != 1
        ax.set_yticks(yt, fill("", length(yt)))
        yl = ax.get_ylim()

    else
        ax.set_ylabel("Depth [m]", fontsize = 15)
        yt = convert(Vector{Int}, yt)
        ax.set_yticks(yt, yt, fontsize = 12)
        yl = ax.get_ylim()
    end
    if i == 11
        tx.set_yticks(depths, string.(keys(locs)), fontsize = 12)
        tx.set_ylabel("Sediment Core", fontsize = 15)
    else
        tx.set_yticks(depths, fill("", length(depths)))
end
    tx.set_ylim(yl)
    title(string(i), fontsize = 15)
    ax.set_xticks([-1, 0, 1], [-1, 0, 1], fontsize = 12)
end
tight_layout()

savefig(plotsdir("U.png"))

#=
# ================= FIGURE X: Comparison with EN4  ============= #
path = "/home/brynn/Code/oceanFP/data/EN4/analyses"
nct = NCDataset(joinpath(path, readdir(path)[1]))
en4lat = nct["lat"][:] #degN
en4lon = nct["lon"][:] #degE
en4depth = nct["depth"][:]
dind = findmin(abs.(en4depth .- 200))[2]
close(nct)
mat = Array{Union{Missing, Float32}}(undef, length(en4lon),length(en4lat), length(readdir(path)), 2)
t = Vector{DateTime}(undef, length(readdir(path)))
for (i,f) in enumerate(readdir(path))
    println("Opening file " * string(i) * " out of " * string(length(readdir(path))))
    nc = NCDataset(joinpath(path, f)) 
    t[i] = nc["time"][1]
    mat[:,:, i, 1] = nc["temperature"][:, :, dind, 1]
    mat[:,:, i, 2] = nc["salinity"][:, :, dind, 1]
    close(nc)
end

tdec = 1900:10:2010
en4 = cat([mean(mat[:, :, (i-1)*12*10+1:i*12*10, :], dims = 3)[:, :, 1,:] for i in 1:length(tdec)]..., dims = 4)
en4[ismissing.(en4)] .= NaN

en4 .-= mean(en4[:, :, :, begin:7], dims = 4)

inds = collect(1900yr:10yr:1970yr)

lev = -1.5:0.2:1.5

proj = ccrs.Orthographic(central_longitude=-80+55,central_latitude = 35)
noproj = ccrs.PlateCarree(central_longitude = 0)
fig, axs = subplots(nrows = length(solutions) + 1, ncols = length(inds),
                    subplot_kw = Dict("projection"=>proj), figsize = (15, 10),
                    constrained_layout = true)
#(80,35,-80,30,5)
lon_box = [-80, 30]
lat_box = [35, 80]
for (i, sol) in enumerate([solutions..., en4])
    Tu = sol isa solution ? Array(sol.ũ.dims[1]) : tdec .* yr
    Ty = sol isa solution ? Array(sol.y.dims[1]) : tdec .* yr 
    anom_index = findall(x->1900yr<x<1970yr, Tu)

    #remove 1850-1970anomaly from value we are plotting 
    if sol isa solution
        θ_ = [ustrip.(value.(x)) for x in sol.θ]
        θ = [x .- (sum(θ_[anom_index]) ./ length(anom_index)) for x in θ_]
    else
        θ = sol
        θ .-= (sum(θ[:, :,1, anom_index], dims = 3) ./ length(anom_index))
    end
            
    for (ii, y) in enumerate(inds)
        
        ind = findmin(abs.(Tu .- y))[2]
        @show ind
        if sol isa solution && y ∈ Ty
            plotme = γreshape(θ[ind], sol.γ)' #100yr at 10yr res 
        elseif sol isa Array
            plotme = θ[:, :,1, ind]
            plotme = plotme'
        end
        if (sol isa solution && y ∈ Ty) || (sol isa Array)
            #ax =  axs[(ii-1) * (length(solutions)+1)+ i]
            ax = axs[i-1, ii-1]
            ax_hdl = ax.plot(orthographic_axes(lat_box..., lon_box..., 5)...,
                             color="black", linewidth=0.5,
                             transform=noproj)
            tx_path = ax_hdl[0]._get_transformed_path()
            path_in_data_coords, _ = tx_path.get_transformed_path_and_affine()
            polygon1s = mpath.Path(path_in_data_coords.vertices)
            ax.set_boundary(polygon1s)
            ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "110m", edgecolor="k", facecolor="gray"))

            if i == 1
                ax.set_title(string(convert(Int64, ustrip(Tu[ind]))), fontsize = 20, fontweight = "bold")
            end
            ax.set_extent([lon_box..., lat_box[1]- 10, lat_box[2]], noproj)
            
            #add a cyclic point
            lon = sol isa solution ? ustrip.(sol.γ.lon) : en4lon
            #plotme, lon_ = cu.add_cyclic_point(plotme, coord = lon)
            lat = sol isa solution ? ustrip.(sol.γ.lat) : en4lat
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
    savefig(plotsdir("surfacesol" * suffix * ".png"), bbox_inches = "tight")
end
=# 
=
