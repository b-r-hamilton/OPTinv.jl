#=
REGION MEAN
=#
if ! @isdefined(solutions)
    include("ex3.transientinversion_abstract.jl")
end

# =================== MEAN OF N. ATL BOX, for each sol'n, compared to LMR and OC2k = #
Tm1 = ustrip.(minimum([sol.y.dims[1][1] for sol in solutions]))
Tm2 = ustrip.(maximum([sol.y.dims[1][end-1] for sol in solutions]))
regionmeanindices = γbox(oldc.γ, 49, 89, 309, 21)
figure()

for sol in solutions
    sol_anom_ind = findall(x->1850yr<x<1970yr, Array(sol.ũ.dims[1]))
    θbox = estimate(sol.ũ, sol.spatialmodes, :θ, spatialinds = regionmeanindices, rolling = 2)
    μ = mean(θbox.v[sol_anom_ind])
    θbox = DimEstimate(θbox.v .- μ, θbox.C, θbox.dims)
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

lmrdataset = loadLMR("sst")
lmrtime = lmrdataset["time"][:]; lmrlon = lmrdataset["lon"][:]; lmrlat = lmrdataset["lat"][:]

for i in 1:20
    lmrsst = makeNaN(lmrdataset["sst"][:, :, i, :])
    anomindex = findall(x->1850 < x < 1970, year.(lmrtime))
    lmrsst = lmrsst .- mean(lmrsst[:, :, anomindex], dims = 3)[:, :, 1]

    lmrtimeroll = rolling(mean, year.(lmrtime), 50)
    lmrglobal = rolling(mean, [NaNMath.mean(lmrsst[:, :, i]) for i in 1:size(lmrsst)[3]], 50)

    lmrbox = rolling(mean, boxmean(lmrsst, lmrlat, lmrlon, 50, 90, 360-50, 20), 50)
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

xlim(Tm1, Tm2)
xticks(800:200:1800, fontsize = 12) 
ylim(-0.5, 0.9)
yticks(-0.5:0.25:0.75, fontsize = 12)
tight_layout()
savefig(plotsdir("meants" * suffix * ".png"))


# ========== NORDIC SEA V. SPNA REGION MEAN  ==================== #

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
    xlim(Tm1, Tm2)
    ylabel("Temperature Anomaly [K]",fontsize = 15)
    if i == 2 xlabel("Time [years CE]", fontsize = 15) end
 
    yticks(-0.5:0.25:0.75, fontsize = 12)
    ylim(-0.5, 0.8)
    xticks(800:200:1800, fontsize = 12)
    text(x = 810, y = -0.48, s = ["A", "B"][i], fontsize = 30, weight = "bold")
    twin2 = ax.twinx()
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
