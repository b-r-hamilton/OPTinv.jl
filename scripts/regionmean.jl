#=
REGION MEAN
=#
if ! @isdefined(solutions)
    include("transientinversion.jl")
end
using UnitfulLinearAlgebra, LinearAlgebra, PaleoData, NaNMath, DateFormats, Dates, RollingFunctions

# =================== MEAN OF N. ATL BOX, for each sol'n, compared to LMR and OC2k = #
Tm1 = ustrip.(minimum([sol.y.dims[1][1] for sol in solutions]))
Tm2 = ustrip.(maximum([sol.y.dims[1][end-1] for sol in solutions]))
regionmeanindices = γbox(solutions[1].γ, 49, 89, 309, 21)
fig, ax1 = subplots(figsize = (8,4))

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
    @show findmax(θbox.x[pre1900inds])[1] - findmin(θbox.x[inds])[1]

    println("LLS, LIA cooling")

    #we need to subset to the time period we want to compute the slope over 
    liainds = findall(x->sol.ũ.dims[1][inds][minv[2]] > x > sol.ũ.dims[1][pre1900inds][max[2]], Array(sol.ũ.dims[1]))
    x = UnitfulMatrix(Array(sol.ũ.dims[1][liainds]))
    C = UnitfulMatrix(parent(θbox.C[liainds, liainds]), unitrange(θbox.C)[liainds], unitdomain(θbox.C)[liainds])
    y = UnitfulMatrix(θbox.v[liainds])
    lls, C = linearleastsquares(x, y, C=C)
    @show lls[1] * 1000 
    @show sqrt.(diag(C))[1] .* 1000

    plot(θbox.x[inds], color = sol.color, label = sol.name, lzorder = 3, fbzorder = 0,lwcentral = 3,lwedges = 1)
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
oc2kzorder = 2
oc2kcolor = "gray"
hlines(y = y_ocean2k, xmin = t_ocean2k .- 100, xmax = t_ocean2k .+ 100, color = oc2kcolor, linewidth = 3, zorder = oc2kzorder)
hlines(y = y_ocean2k .+ ystd_ocean2k[2, :], xmin = t_ocean2k .- 100, xmax = t_ocean2k .+ 100, color = oc2kcolor,zorder = oc2kzorder)
hlines(y = y_ocean2k .- ystd_ocean2k[1, :], xmin = t_ocean2k .- 100, xmax = t_ocean2k .+ 100, color = oc2kcolor,zorder = oc2kzorder)
vlines(x = t_ocean2k .- 100, ymin = y_ocean2k .- ystd_ocean2k[1, :], ymax = y_ocean2k .+ ystd_ocean2k[2,:], color = oc2kcolor,zorder = oc2kzorder)
vlines(x = t_ocean2k .+ 100, ymin = y_ocean2k .- ystd_ocean2k[1, :], ymax = y_ocean2k .+ ystd_ocean2k[2,:], color = oc2kcolor,zorder = oc2kzorder)
[fill_between(x = [t-100, t+100], y1 = y .- ystd1, y2 = y .+ ystd2, color = oc2kcolor,alpha = 0.5, zorder = oc2kzorder) for (t, y, ystd1, ystd2) in zip(t_ocean2k, y_ocean2k, ystd_ocean2k[1, :], ystd_ocean2k[2, :])]

lmrdataset = loadLMR("sst")
lmrtime = lmrdataset["time"][:]; lmrlon = lmrdataset["lon"][:]; lmrlat = lmrdataset["lat"][:]

month.(lmrtime)

for i in 1:20
    lmrsst = makeNaN(lmrdataset["sst"][:, :, i, :])
    anomindex = findall(x->1850 < x < 1970, year.(lmrtime))
    lmrsst = lmrsst .- mean(lmrsst[:, :, anomindex], dims = 3)[:, :, 1]

    lmrtimeroll = rolling(mean, year.(lmrtime), 50)
    lmrglobal = rolling(mean, [NaNMath.mean(lmrsst[:, :, i]) for i in 1:size(lmrsst)[3]], 50)

    lmrbox = rolling(mean, boxmean(lmrsst, lmrlat, lmrlon, 50, 90, 360-50, 20), 50)
    plot(DimArray((lmrbox  .- mean(lmrbox))* K, Ti(lmrtimeroll)), color = "black", label = "LMR", zorder = 1, linewidth = 0.5, linestyle = "solid")
end

ylabel("Surface Temperature Anomaly [K]", fontsize = 15)
xlabel("Time [yr CE]", fontsize = 15)

xlim(Tm1, Tm2)
xticks(800:200:1800, fontsize = 12) 
ax1.set_ylim(-0.42, 0.75)
yticks(-0.25:0.25:0.75, fontsize = 12)
tight_layout()

inset_ax = ax1.inset_axes([1550,0.4,300,0.3], transform = ax1.transData ) 
#subplot
hadisst = loadHadISST()
hadsst = makeNaN(hadisst["tos"][:, :, :])
hadlat = hadisst["latitude"][:]; hadlon = hadisst["longitude"][:]; hadtime = hadisst["time"][:]
hadroll(x) = rolling(mean, x, 50 * 12)
hadisst_gm = hadroll([NaNMath.mean(hadsst[:, :, i]) for i in 1:size(hadsst)[3]])
hadisst_bm = hadroll(boxmean(hadsst, hadlat, hadlon, 50, 90, -50, 20))
hadtimeroll = hadroll(yeardecimal.(hadtime))
# hadisst_bm = boxmean(hadsst, hadlat, hadlon, 50, 90, -50, 20)
# hadtimeroll = yeardecimal.(hadtime)
hadtimeanom = findall(x->1850<x<1970, hadtimeroll)
inset_ax.plot(hadtimeroll, hadisst_bm .- mean(hadisst_bm[hadtimeanom]), color = "black", label = "HadISST", linewidth = 3, linestyle = "dashed")
for sol in solutions
    sol_anom_ind = findall(x->1850yr<x<1970yr, Array(sol.ũ.dims[1]))
    θbox = estimate(sol.ũ, sol.spatialmodes, :θ, spatialinds = regionmeanindices, rolling = 2)
    inds = findall(x->x∈sol.y.dims[1], Array(sol.ũ.dims[1]))
    μ = mean(θbox.v[sol_anom_ind])
    θbox = DimEstimate(θbox.v .- μ, θbox.C, θbox.dims)
    t = ustrip.(Array(θbox.x[inds].dims[1]))
    y = ustrip.(θbox.v[inds])
    yunc = ustrip.(Array(uncertainty.(θbox.x[inds])))
    inset_ax.plot(t, y, color = sol.color, label = sol.name, zorder = 10000)
    inset_ax.fill_between(x = t, y1 = y .- yunc, y2 = y .+ yunc, zorder = 10000, color = sol.color, alpha = 0.3)
end

inset_ax.set_xlim(1875,1970)
inset_ax.set_ylim(-0.4, 0.5)
inset_ax.set_ylabel("T [K]")

xl = inset_ax.get_xlim()
yl = inset_ax.get_ylim()

ax1.hlines(xmin = xl[0], xmax = xl[1], y = yl, color = "black")
ax1.vlines(ymin = yl[0], ymax = yl[1], x = xl, color = "black")
savefig(plotsdir("meants" * suffix * ".png"), dpi = 600)

# ================== LABRADOR SEA =============================== #
figure(figsize = (8, 8))
lab_lon = (360-60, 360-45)
lab_lat = (55, 64)
df = CSV.read(OPTinv.datadir("TLS1998Fig7Digitized.csv"), DataFrame)
df = df[sortperm(df[!, "Year"]), :]

tls_anom_ind = findall(x->x == 70, df[!, "Year"])
tls_sal = df[!, "Salinity"] .- df[tls_anom_ind, "Salinity"]
tls_temp = df[!, "Temperature"] .- df[tls_anom_ind, "Temperature"]
tls_time = 1900 .+ df[!, "Year"]

labmeanindices = γbox(solutions[1].γ, 55, 64, 360-60, 360-45)
subplot(2,1,1)
text(x = 1901, y = 0.75, s = "A", fontsize = 30, weight = "bold")
grid()
for sol in solutions
    sol_anom_ind = findall(x->x==1970yr, Array(sol.ũ.dims[1]))
    labθ = estimate(sol.ũ, sol.spatialmodes, :θ, spatialinds = labmeanindices, rolling = 2)
    labθ = DimEstimate(labθ.v .- value(labθ.x[At(1970yr), :]), labθ.C, labθ.dims)
    plot(labθ.x, color = sol.color, label = sol.name, lzorder = 2)
end
    
plot(1900 .+ df[!, "Year"], tls_temp, color = "black", label = "TLS1998", zorder = 1, linestyle = "solid", linewidth = 3)
xlim(1900, 1970)
ylim(-1,1)
xticks(1900:10:1970, fontsize = 12)
yticks(-1:0.5:1, fontsize = 12)
ylabel("Temperature Anomaly [K]", fontsize = 15)

subplot(2,1,2)
text(x = 1901, y = 0.045, s = "B", fontsize = 30, weight = "bold")
for sol in solutions
    sol_anom_ind = findall(x->x==1970yr, Array(sol.ũ.dims[1]))
    labS = estimate(sol.ũ, sol.spatialmodes, :δ, spatialinds = labmeanindices, rolling = 2) 
    labS = DimEstimate((labS.v .- value(labS.x[At(1970yr), :])) ./ 0.5255, labS.C, labS.dims)
    plot(labS.x, color = sol.color, label = sol.name, lzorder = 2)
end

plot(1900 .+ df[!, "Year"], tls_sal, color = "black", label = "TLS1998", linestyle = "solid", linewidth = 3)

xlim(1900, 1970)
grid()
ylim(-0.06, 0.06)
ylabel("Salinity Anomaly [g/kg]", fontsize = 15)
xlabel("Time [years CE]", fontsize = 15)
xticks(1900:10:1970, fontsize = 12)
yticks(-0.06:0.03:0.06, fontsize = 12)
tight_layout()
savefig(plotsdir("labcomp" * suffix * ".png")) 

#=
# ================== OPT-3, OPT-11 vs. GH19 ================== #
Tm1 = ustrip.(minimum([sol.y.dims[1][1] for sol in solutions]))
Tm2 = ustrip.(maximum([sol.y.dims[1][end-1] for sol in solutions]))
regionmeanindices = γbox(solutions[1].γ, 49, 89, 309, 21)
fig, ax1 = subplots(figsize = (8,4))

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
    @show findmax(θbox.x[pre1900inds])[1] - findmin(θbox.x[inds])[1]

    println("LLS, LIA cooling")
    #we need to subset to the time period we want to compute the slope over 
    liainds = findall(x->sol.ũ.dims[1][inds][minv[2]] > x > sol.ũ.dims[1][pre1900inds][max[2]], Array(sol.ũ.dims[1]))
    x = UnitfulMatrix(Array(sol.ũ.dims[1][liainds]))
    C = UnitfulMatrix(parent(θbox.C[liainds, liainds]), unitrange(θbox.C)[liainds], unitdomain(θbox.C)[liainds])
    y = UnitfulMatrix(θbox.v[liainds])
    lls, C = linearleastsquares(x, y, C=C)
    @show lls[1] * 1000 
    @show sqrt.(diag(C))[1] .* 1000

    plot(θbox.x[inds], color = sol.color, label = sol.name, lzorder = 3, fbzorder = 0,lwcentral = 3,lwedges = 1)
    newcolor = sol.color == "red" ? "pink" : "aqua"
    println()
end

using NCDatasets, GH19 
exp = explist()[end]
filename = download_exp(exp, true)
tmi = NCDataset(TMI.download_ncfile(TMIversion()))
lat = tmi["lat"][:]
lon = tmi["lon"][:]
natl = tmi["d_NATL"][:, :]
natl[natl .== 0] .= NaN
ant =  tmi["d_ANT"][:, :]
ant[ant .== 0] .= NaN
glob = tmi["d_GLOBAL"][:,:]
glob[glob .== 0] .= NaN               
nc = NCDataset(filename)
θGH19 = nc["theta"][:, 1, :, :]
tGH19 = nc["year"][:]
GHanomindex = findall(x->1850 < x < 1970, tGH19)
θGH19NATL = [NaNMath.mean(θGH19[i,:, :] .* natl') for i in 1:size(θGH19)[1]]
θGH19ANT = [NaNMath.mean(θGH19[i,:, :] .* ant') for i in 1:size(θGH19)[1]]
θGH19GLOB = [NaNMath.mean(θGH19[i,:, :] .* glob') for i in 1:size(θGH19)[1]]
ind1970 = findall(x->x==1970, tGH19)[1]
ind1880 = findall(x->x==1880, tGH19)[1]
ind1780 = findall(x->x==1780, tGH19)[1]    
θGH19NATL[ind1970]-θGH19NATL[ind1780]
θGH19NATL[ind1970]-θGH19NATL[ind1780]


plot(tGH19, θGH19NATL .- mean(θGH19NATL[GHanomindex]), color = "green", linestyle = "solid", linewidth = 5)
plot(tGH19, θGH19GLOB .- mean(θGH19GLOB[GHanomindex]), color = "black", linestyle = "solid", linewidth = 5)



tGH19
θGH19NATL
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
oc2kzorder = 2
oc2kcolor = "gray"
hlines(y = y_ocean2k, xmin = t_ocean2k .- 100, xmax = t_ocean2k .+ 100, color = oc2kcolor, linewidth = 3, zorder = oc2kzorder)
hlines(y = y_ocean2k .+ ystd_ocean2k[2, :], xmin = t_ocean2k .- 100, xmax = t_ocean2k .+ 100, color = oc2kcolor,zorder = oc2kzorder)
hlines(y = y_ocean2k .- ystd_ocean2k[1, :], xmin = t_ocean2k .- 100, xmax = t_ocean2k .+ 100, color = oc2kcolor,zorder = oc2kzorder)
vlines(x = t_ocean2k .- 100, ymin = y_ocean2k .- ystd_ocean2k[1, :], ymax = y_ocean2k .+ ystd_ocean2k[2,:], color = oc2kcolor,zorder = oc2kzorder)
vlines(x = t_ocean2k .+ 100, ymin = y_ocean2k .- ystd_ocean2k[1, :], ymax = y_ocean2k .+ ystd_ocean2k[2,:], color = oc2kcolor,zorder = oc2kzorder)
[fill_between(x = [t-100, t+100], y1 = y .- ystd1, y2 = y .+ ystd2, color = oc2kcolor,alpha = 0.5, zorder = oc2kzorder) for (t, y, ystd1, ystd2) in zip(t_ocean2k, y_ocean2k, ystd_ocean2k[1, :], ystd_ocean2k[2, :])]

ylabel("Surface Temperature Anomaly [K]", fontsize = 15)
xlabel("Time [yr CE]", fontsize = 15)

xlim(Tm1, Tm2)
xticks(800:200:1800, fontsize = 12) 
ax1.set_ylim(-0.7, 1.1)
yticks(-0.5:0.25:0.75, fontsize = 12)
tight_layout()
savefig(plotsdir("meants_gh19" * suffix * ".png"), dpi = 600)
=#

