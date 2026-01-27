#=
REGION MEAN
=#

if ! @isdefined(solutions)
    include("transientinversion.jl")
end
using UnitfulLinearAlgebra, LinearAlgebra, PaleoData, NaNMath, DateFormats, Dates, RollingFunctions, SkipNan, PythonCall, CSV, DataFrames

# =================== MEAN OF N. ATL BOX, for each sol'n, compared to LMR and OC2k = #
Tm1 = ustrip.(minimum([sol.y.dims[1][1] for sol in solutions]))
Tm2 = ustrip.(maximum([sol.y.dims[1][end-1] for sol in solutions]))
ΔT = unique(diff(Array(solutions[1].y.dims[1])))[1] #ASSUMES EACH SOL HAS SAME TIMESTEP
regionmeanindices = γbox(solutions[1].γ, 49, 89, 309, 21)
ca = vec(cellarea(solutions[1].γ).tracer)[regionmeanindices]
∑ca = sum(ca)
ca ./= sum(ca)
fig, ax1 = subplots(figsize = (8,4))
roll = 2
rollyrs = isnothing(roll) ? ΔT : ((2 * roll) + 1) * ΔT
rollyrs = convert(Int64, ustrip(rollyrs))
for sol in solutions
    sol_anom_ind = findall(x->1850yr<x<1970yr, Array(sol.ũ.dims[1]))
    θbox = estimate(sol.ũ, sol.spatialmodes, :θ, spatialinds = regionmeanindices, rolling = roll, weights = ca)
    μ = mean(θbox.v[sol_anom_ind])
    θbox = DimEstimate(θbox.v .- μ, θbox.C, θbox.dims)
    inds = findall(x->x∈sol.y.dims[1], Array(sol.ũ.dims[1]))
    pre1900inds = intersect(inds, findall(x->1000yr<x<1900yr, Array(sol.ũ.dims[1])))
    minv = findmin(θbox.x[inds])
    maxv = findmax(θbox.x[pre1900inds])
    println(sol.name)
    println("Time of minima: " * string(sol.ũ.dims[1][inds][minv[2]])) 
    println("Time of maxima (pre1900): "*string(sol.ũ.dims[1][pre1900inds][maxv[2]]))
    
    #we need to subset to the time period we want to compute the slope over 
    liainds = findall(x->sol.ũ.dims[1][inds][minv[2]] ≥ x ≥ sol.ũ.dims[1][pre1900inds][maxv[2]], Array(sol.ũ.dims[1]))
    x = UnitfulMatrix(Array(sol.ũ.dims[1][liainds]))
    C = UnitfulMatrix(parent(θbox.C[liainds, liainds]), unitrange(θbox.C)[liainds], unitdomain(θbox.C)[liainds])
    y = UnitfulMatrix(θbox.v[liainds])
    lls, C = linearleastsquares(x, y, C=C)
    println("LIA cooling rate: " * string(round(lls[1] * 1000, sigdigits = 2)) * "±" * string(round(ustrip(sqrt.(diag(C))[1]) .* 1000, sigdigits = 2)) * " K/kyr")
    println()
            
    plot(θbox.x[inds], color = sol.color, label = sol.name, lzorder = 3, fbzorder = 0,lwcentral = 3,lwedges = 1)
    newcolor = sol.color == "red" ? "pink" : "aqua"

end
oc2k_binned, oc2k_N, oc2k_ages, oc2k_binnedv, oc2k = loadOcean2kBinned()

#=
f1(x,f) = [f(skipnan(x[i, :])) for i in 1:size(x)[1]]
f2(x,f) = [f(skipnan(x[:, i])) for i in 1:size(x)[2]]



oc2k_anomrem = oc2k_binned .- repeat(f2(oc2k_binned, mean)', inner = (10, 1))
globalanom = f1(weighted_regional_mean, sum)

figure();plot(t_ocean2k, globalanom);ylim(-3,3)
[convert(Vector{Int64}, vec(oc2k[r .* "d"])) for r in regions]

tlls =  convert(Vector{Float64}, collect(t_ocean2k)) ./ 1000
weighted_regional_mean

function lls_remnan(t,y)
    nanind = findall(x->!isnan(x),y)
    return linearleastsquares(t[nanind], y[nanind])
end
for i in 1:57
    @show i 
    lls_remnan(tlls, oc2k_anomrem[:,i])
end

=#

t_ocean2k = collect(100:200:1900)
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

if ! @isdefined(lmrdataset)
    lmrdataset = loadLMR("sst")
    lmrtime = lmrdataset["time"][:]; lmrlon = lmrdataset["lon"][:]; lmrlat = lmrdataset["lat"][:]
end

##= SUBMISSION 1 LMR PLOTTING 
for i in 1:20
    lmrsst = makeNaN(lmrdataset["sst"][:, :, i, :])
    anomindex = findall(x->1850 < x < 1970, year.(lmrtime))
    lmrsst = lmrsst .- mean(lmrsst[:, :, anomindex], dims = 3)[:, :, 1]


    lmrtimeroll = rolling(mean, year.(lmrtime), rollyrs)
    lmrglobal = rolling(mean,
                        [NaNMath.mean(lmrsst[:, :, i]) for i in 1:size(lmrsst)[3]],
                        rollyrs)
    lmrbox = rolling(mean, boxmean(lmrsst, lmrlat, lmrlon, 50, 90, 360-50, 20),
                     rollyrs)
    plot(DimArray((lmrbox  .- mean(lmrbox))* K, Ti(lmrtimeroll)), color = "black", label = "LMR", zorder = 1, linewidth = 0.5, linestyle = "solid")
end


# = SUBMISSION 2 LMR PLOTTING
wet = (!).(ismissing.(lmrdataset["sst"][:, :, 1, 1]))
I = cartesianindex(wet)
area = computearea(lmrlon, lmrlat, I)
natlbox =  vcat([[CartesianIndex(i, j) for i in findall(x->20>x|| x>360-50, lmrlon)] for j in findall(x->50<x<90, lmrlat)]...)
notnatlbox = [i for i in I if i ∉ natlbox]
area[notnatlbox] .= 0 

areaweight = area ./ sum(area) .* wet
lmraw = lmrdataset["sst"][:, :, :, :] .* areaweight
lmrsstμ = mean(lmraw, dims = 3)[:, :, 1, :]
lmrsststd = std(lmraw, dims = 3)[:, :, 1, :]

sst_aw_mean = rolling(mean, [sum(skipmissing(lmrsstμ[:, :, i])) for i in 1:2001], 50)
sst_aw_std =  rolling(mean, [sum(skipmissing(lmrsststd[:, :, i])) for i in 1:2001], 50)
rlmrtime = rolling(mean, year.(lmrtime), 50)
lmr_anom_ind = findall(x->1850<x<1970, rlmrtime)
sst_aw_mean .-= mean(sst_aw_mean[lmr_anom_ind])

#plot(rlmrtime, sst_aw_mean, color = "black")
#fill_between(rlmrtime, y1 = sst_aw_mean .- sst_aw_std, y2 = sst_aw_mean .+ sst_aw_std, color = "black", alpha = 0.2)


 tind = findall(x->800<x<1800, year.(lmrtime)) #ocean2k indices 
 lls, llsC = linearleastsquares(convert(Vector{Float64}, year.(lmrtime))[tind],
                               sst_aw_mean[tind], C = diagm(sst_aw_std[tind].^2))

#plot(year.(lmrtime)[tind], year.(lmrtime)[tind] .* lls[1] .+ lls[2], color = "red")
#ylim(-0.4, 0.75)
println("LMR cooling rate: " *string(round(lls[1,1] * 1000, sigdigits = 2))  * string("±") * string(round(sqrt(llsC[1,1]) * 1000, sigdigits = 2)) * " °C/kyr")


ylabel("Surface Temperature Anomaly [K]", fontsize = 15)
xlabel("Time [yr CE]", fontsize = 15)

xlim(Tm1, Tm2)
xticks(800:200:1800, fontsize = 12) 
ax1.set_ylim(-0.42, 0.75)
yticks(-0.25:0.25:0.75, fontsize = 12)
tight_layout()

inset_ax = ax1.inset_axes([1550,0.4,300,0.3], transform = ax1.transData ) 
#subplot
if ! @isdefined(hadisst)
    hadisst = loadHadISST()
    hadsst = makeNaN(hadisst["tos"][:, :, :])
    hadlat = hadisst["latitude"][:]; hadlon = hadisst["longitude"][:]; hadtime = hadisst["time"][:]
end


if !isnothing(rolling)
    hadroll(x) = rolling(mean, x, 50 * 12)
else
    hadroll(x) = rolling(mean, x, ustrip(rollyrs) * 12)
end

hadisst_gm = hadroll([NaNMath.mean(hadsst[:, :, i]) for i in 1:size(hadsst)[3]])
hadisst_bm = hadroll(boxmean(hadsst, hadlat, hadlon, 50, 90, -50, 20))
hadtimeroll = hadroll(yeardecimal.(hadtime))
# hadisst_bm = boxmean(hadsst, hadlat, hadlon, 50, 90, -50, 20)
# hadtimeroll = yeardecimal.(hadtime)
hadtimeanom = findall(x->1850<x<1970, hadtimeroll)
inset_ax.plot(hadtimeroll, hadisst_bm .- mean(hadisst_bm[hadtimeanom]), color = "black", label = "HadISST", linewidth = 3, linestyle = "dashed")
for sol in solutions
    sol_anom_ind = findall(x->1850yr<x<1970yr, Array(sol.ũ.dims[1]))
    θbox = estimate(sol.ũ, sol.spatialmodes, :θ, spatialinds = regionmeanindices, rolling = roll)
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
    labθ = estimate(sol.ũ, sol.spatialmodes, :θ, spatialinds = labmeanindices, rolling = roll)
    labθ = DimEstimate(labθ.v .- value(labθ.x[Near(1970yr), :]), labθ.C, labθ.dims)
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
    labS = estimate(sol.ũ, sol.spatialmodes, :δ, spatialinds = labmeanindices, rolling = roll) 
    labS = DimEstimate((labS.v .- value(labS.x[Near(1970yr), :])) ./ 0.5255, labS.C, labS.dims)
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

# SUBMISSION 2, d18Osw results
fig, ax1 = subplots(figsize = (8,4))
for sol in solutions
    sol_anom_ind = findall(x->1850yr<x<1970yr, Array(sol.ũ.dims[1]))
    δbox = estimate(sol.ũ, sol.spatialmodes, :δ, spatialinds = regionmeanindices, rolling = roll, weights = ca)

    μ = mean(δbox.v[sol_anom_ind])
    δbox = DimEstimate(δbox.v .- μ, δbox.C, δbox.dims)
    inds = findall(x->x∈sol.y.dims[1], Array(sol.ũ.dims[1]))
    #Sbox = 2.796 .* δbox.x .+ 34.38 
    plot(δbox.x[inds], color = sol.color, label = sol.name, lzorder = 3, fbzorder = 0,lwcentral = 3,lwedges = 1)
end
yl = ax1.get_ylim()
yl = (yl[0], yl[1])
ylnew = yl .* 2.796 #.+ 34.38
ax2 = ax1.twinx()
ax2.set_ylim(ylnew)
ax1.set_ylabel("Surface δ¹⁸O" * L"_{\mathrm{seawater}}" * " Anomaly [‰]", fontsize = 15)
ax1.set_xlabel("Time [yr CE]", fontsize = 15)
ax2.set_ylabel("Salinity Anomaly [psu]", fontsize = 15)
xlim(Tm1, Tm2)
ax1.set_xticklabels(800:200:1800, fontsize = 12)
ax1.set_yticklabels(-0.03:0.01:0.04,fontsize = 12)
ax2.set_yticklabels(-0.075:0.025:0.125,fontsize = 12) 
tight_layout()
savefig(plotsdir("meantsd18O" * suffix * ".png"), dpi = 600)

# ============== WHAT DOES RECONSTRUCTION LOOK LIKE WITH JUST MODE 1 ============== #

#find indices associated with NOT mode 1 
ũcp = deepcopy(oldc.ũ.x)
#zero those values out 
ũcp[:, At(2:11), :] .*= 0
zind = vec(ũcp .== 0)
#make adjusted covariance matrix where we zero out the modes not associated with Mode 1 
Ccp = deepcopy(oldc.ũ.C);
Ccp[zind, zind] .*= 0;

sol_anom_ind = findall(x->1850yr<x<1970yr, Array(oldc.ũ.dims[1]))
inds = findall(x->x∈oldc.y.dims[1], Array(oldc.ũ.dims[1]))
ũmode1 = DimEstimate(UnitfulMatrix(vec(value.(ũcp))), Ccp, oldc.ũ.dims)
θboxmode1 = estimate(ũmode1, oldc.spatialmodes, :θ, spatialinds = regionmeanindices, rolling = roll, weights = ca)
μ = mean(θboxmode1.v[sol_anom_ind])
θboxmode1 = DimEstimate(θboxmode1.v .- μ, θboxmode1.C, θboxmode1.dims)

θboxold = estimate(oldc.ũ, oldc.spatialmodes, :θ, spatialinds = regionmeanindices, rolling = roll, weights = ca)
μ = mean(θboxold.v[sol_anom_ind])
θboxold = DimEstimate(θboxold.v .- μ, θboxold.C, θboxmode1.dims)

figure()
plot(θboxmode1.x[inds], color = "red", label = "Mode 1", lzorder = 3, fbzorder = 0,lwcentral = 3,lwedges = 1)
plot(θboxold.x[inds], color = "blue", label = "old", lzorder = 3, fbzorder = 0,lwcentral = 3,lwedges = 1)

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

