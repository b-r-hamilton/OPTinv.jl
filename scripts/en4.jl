
#=
Right now, this code can compute T at bottom of mixed layer for EN4 profiles
according to the Holte and Talley algorithm

How do we compare this to my reconstructions?

Data Coverage is very poor pre-1940, and I only can compute MLD T for good data

Here's an idea
1. Compute March bottom of MLD T from mean EN4 gridded product
2. For every available record, compute its bottom of MLD T and subtract it from nearest EN4 gridded product
3. Compute your anomaly maps from the EN4 time period (1900-pres)
4. Overlay with March bottom of MLD T and S anom. obs. as scatter points 
=#

import Pkg
Pkg.activate("../")
using NCDatasets, Dates, CSV, DataFrames, Statistics, OPTinv, GibbsSeaWater, PyCall, PyPlot, NaNMath, DimensionalData, DrWatson, JLD2
#ran Conda.pip_interop(true); Conda.pip("install", "holteandtalley")
hat = pyimport("holteandtalley")
cm = pyimport("cmocean.cm") 
ccrs = pyimport("cartopy.crs")

# ============== 1.  Compute climatological mean bottom of MLD T and S ========== #
savename = datadir("EN4_march_mld.jld2")
    lonbox = [-80, 30]
    latbox = [35, 80]

if !isfile(savename)
    println("Computing ML T and S, will take approx. 30 mins.") 
    path = "/home/brynn/Code/oceanFP/data/EN4/analyses"
    
    files = readdir(path)
    nc = NCDataset(joinpath(path, files[1]))

    lat = nc["lat"][:]
    lon = nc["lon"][:]
    z = nc["depth"][:]
    p = hcat([gsw_p_from_z.(-1 * z, l) for l in lat]...)' #lat x depth

    pstring(s::String) = (parse(Int, s[end-8:end-5]), parse(Int, s[end-4:end-3]));pstring(files[1])
    ft = pstring.(files)
    ft_times = [Dates.Date(f[1], f[2], 1) for f in ft]
    march = findall(x->x == 3, month.(ft_times))

    lon_indices = findall(x->360-80 < x || x < 30, lon)
    lat_indices = findall(x->latbox[1] < x < latbox[2], lat)
    p = p[lat_indices, :]

    mldtemp = Array{Float64}(undef, length(lon_indices), length(lat_indices), length(march)) .* NaN
    mlddepthtemp = Array{Float64}(undef, length(lon_indices), length(lat_indices), length(march)) .* NaN
    mldsalinity = Array{Float64}(undef, length(lon_indices), length(lat_indices), length(march)) .* NaN
    mlddepthsalinity = Array{Float64}(undef, length(lon_indices), length(lat_indices), length(march)) .* NaN

    close(nc) 

    #15 seconds per year, for 1900-2022, we have 123 years of data
    #should take 30 minutes to compute 
    for (arrayind, fileind) in enumerate(march)
        nc = NCDataset(joinpath(path, files[fileind]))
        salinity = nc["salinity"][lon_indices, lat_indices, :, 1]
        salinity[ismissing.(salinity)] .= NaN
        salinity = convert(Array{Float64}, salinity)
        temperature = nc["temperature"][lon_indices, lat_indices, :, 1]
        temperature[ismissing.(temperature)] .= NaN
        temperature = convert(Array{Float64}, temperature)
        SA = cat([hcat([gsw_sa_from_sp.(salinity[lon_, lat_, :], p[lat_, :], lon_, lat_) for lat_ in 1:length(lat_indices)]...) for lon_ in 1:length(lon_indices)]..., dims = 3); #0.1 secs
        SA = permutedims(SA, (3,2,1))
        CT = cat([hcat([gsw_ct_from_t.(SA[lon_, lat_, :], temperature[lon_, lat_, :] .- 273.15, p[lat_, :]) for lat_ in 1:length(lat_indices)]...) for lon_ in 1:length(lon_indices)]..., dims = 3);
        CT = permutedims(CT, (3,2,1))
        ρ = cat([hcat([gsw_rho.(SA[lon_, lat_, :], CT[lon_, lat_, :], p[lat_, :]) for lat_ in 1:length(lat_indices)]...) for lon_ in 1:length(lon_indices)]..., dims = 3)
        ρ = permutedims(ρ, (3,2,1))
        nind = (!).(isnan.(temperature))
        
        for lon_ in 1:length(lon_indices)
            for lat_ in 1:length(lat_indices)
                if sum(nind[lon_, lat_, :]) > 3
                    ml = hat.HolteAndTalley(p[lat_, :][nind[lon_, lat_, :]],temperature[lon_, lat_,  :][nind[lon_, lat_, :]] .- 273.15, salinity[lon_, lat_, :][nind[lon_, lat_, :]], ρ[lon_, lat_, :][nind[lon_, lat_, :]])
                    idepthtemp = findmin(abs.(z .- ml.tempMLD))[2]
                    idepthsal = findmin(abs.(z .- ml.salinityMLD))[2]
                    mldtemp[lon_, lat_, arrayind] = temperature[lon_, lat_, idepthtemp]
                    mlddepthtemp[lon_, lat_, arrayind] = z[idepthtemp]
                    mldsalinity[lon_, lat_, arrayind] = salinity[lon_, lat_, idepthsal]
                    mlddepthsalinity[lon_, lat_, arrayind] = z[idepthsal]
                end
            end
        end
        println(string(ft[fileind]) * " processed")
    end

    breakpoint = findmax(diff(lon[lon_indices]))[2]
    lon_rs = vcat(-360 .+ lon[lon_indices][breakpoint+1:end], lon[lon_indices][begin:breakpoint])

    swap(x::Array, breakpoint::Number) = vcat(x[breakpoint+1:end, :, :], x[begin:breakpoint, :, :])
    mldtemp = swap(mldtemp, breakpoint)
    mldsalinity = swap(mldsalinity, breakpoint)
    mlddepthtemp = swap(mlddepthtemp, breakpoint)
    mlddepthsalinity = swap(mlddepthsalinity, breakpoint)
    lat_rs = lat[lat_indices]
    times = ft_times[march]
    jldsave(savename; mldtemp, mldsalinity, mlddepthtemp, mlddepthsalinity, lon_rs, lat_rs, z, times)
else
    println("Loading in ML T and S from file") 
    jld = jldopen(savename)
    mldtemp = jld["mldtemp"]
    mldsalinity = jld["mldsalinity"]
    mlddepthtemp = jld["mlddepthtemp"]
    mlddepthsalinity = jld["mlddepthsalinity"]
    lon_rs = jld["lon_rs"]
    lat_rs = jld["lat_rs"]
    z = jld["z"]
    close(jld)
end

#march climatological MLD T 
mldT = mean(mldtemp, dims = 3)[:, :, 1]
mldS = mean(mldsalinity, dims = 3)[:, :, 1]

# =========== STEP 2: GET MLD T and S from PROFILES ======= #

path = "/home/brynn/Code/oceanFP/data/EN4/profiles"
files = readdir(path)
pstring(s::String) = (parse(Int, s[end-8:end-5]), parse(Int, s[end-4:end-3]));pstring(files[1])
ft = pstring.(files)
ft_times = [Dates.Date(f[1], f[2], 1) for f in ft]
march = findall(x->x == 3, month.(ft_times))

nc = NCDataset(joinpath(path, files[1]))

function record(lat, lon, depth, t, T, S)
    nandices = (!).(ismissing.(depth))
    depth = convert(Vector{Float64}, depth[nandices])
    T[ismissing.(T)] .= NaN
    S[ismissing.(S)] .= NaN
    T = convert(Vector{Float64}, T[nandices]) 
    S = convert(Vector{Float64}, S[nandices])
    p = gsw_p_from_z.(-1*depth, lat)
    SA = gsw_sa_from_sp.(S, p, lon, lat)
    CT = gsw_ct_from_t.(SA, T, p)
    ρ = gsw_rho.(SA, CT, p)
    mls = NaN
    mlt = NaN 
    try
        ml = hat.HolteAndTalley(p, T, S, ρ)
        #find depth index closest to ML 
        mlt = T[findclosest(ml.tempMLD, depth)]
        mls = S[findclosest(ml.salinityMLD, depth)] 
    catch
    end
    return record(lat, lon, depth, t, T, S, p, CT, SA, ρ, mlt, mls)
end

function getrecords(nc::NCDataset, boxes)
    println("accessing record")
    T = nc["TEMP"][:, :]
    lat = nc["LATITUDE"][:]
    lon = nc["LONGITUDE"][:]
    depth = nc["DEPH_CORRECTED"][:, :]
    S = nc["PSAL_CORRECTED"][:,:]
    t = nc["JULD"][:]
    
    ind = NamedTuple{keys(boxes)}((indices(boxes[k]..., lon, lat) for k in keys(boxes)))
    return NamedTuple{keys(boxes)}(([record(lat[i], lon[i], depth[:, i], t[i], T[:, i], S[:, i]) for i in ind[k]] for k in keys(boxes)))
end

boxes = NamedTuple{(:NATL, )}([[lonbox, latbox]])
indices(lonbox, latbox, lon, lat) = intersect(findall(x->lonbox[1]<x<lonbox[2], lon), findall(x->latbox[1]<x<latbox[2], lat))
recs = [getrecords(NCDataset(joinpath(path, f)), boxes) for f in files[march]]
recstime = ft_times[march]

jldsave(datadir("en4_recs.jld2"); recs, recstime)

tstarts = 1899:10:1969
tends = 1910:10:1980
figure();
for (i, (ts, te)) in enumerate(zip(tstarts, tends))
    tind = findall(x->ts < x < te, year.(ft_times[march]))
    recs_of_interest = vcat([r.NATL for r in recs[tind]]...)
    @show length(recs_of_interest)

    lons = []
    lats = []
    anomT = []
    anomS = []
    for r in recs_of_interest
        if !isnan(r.MLDT)
            lonind = findclosest(r.lon, lon_rs)
            latind = findclosest(r.lat, lat_rs)
            push!(lons, r.lon)
            push!(lats, r.lat)
            push!(anomT, r.MLDT +273.15 - mldT[lonind, latind])
            push!(anomS, r.MLDS  - mldS[lonind, latind])
            
        end
    end

ax = subplot(4,2,i;projection = ccrs.PlateCarree())
s = scatter(lons, lats, c = anomT, cmap = cm.balance, vmin = -5, vmax = 5)
#colorbar()
    ax.coastlines()
    title(string(ts) * "-" * string(te))
end
tight_layout()












#=
function compute_ML(T::Vector, S::Vector, lon::Number, lat::Number)
    nandices = (!).(ismissing.(T))
    p_ = vec(p[nandices, At(lat)])
    SA = gsw_sa_from_sp.(S[nandices], p_, lon, lat)
    ρ = gsw_rho.(SA, gsw_ct_from_t.(SA, T[nandices], p_), p_)
    ml = hat.HolteAndTalley(p_, T[nandices], S[nandices], ρ)
    return ml.tempMLD, ml.salinityMLD
end

@time ml = compute_ML(nc["temperature"][lonindex, latindex, :, 1], nc["salinity"][lonindex, latindex, :, 1], nc["lon"][lonindex], nc["lat"][latindex])


figure(); plot(nc["salinity"][lonindex, latindex, :, 1], z); ylim(100, 0)

typeof(nc["temperature"][lonindex, latindex, :, 1])
typeof(z)

for (t, f) in zip(ft_times, files[1])
    SA
    [
end


@time r = record(nc["lat"][latindex], nc["lon"][lonindex], nc["depth"][:], nc["time"][1], nc["temperature"][lonindex, latindex, :, 1], nc["salinity"][lonindex, latindex, :, 1]);


boxes = (
SG = [[-70, -35], [50, 65]], 
SI = [[-30, 0], [50, 65]], 
EG = [[-30, 12], [70, 80]],
)






recs = [getrecords(NCDataset(joinpath(path, f)), boxes) for f in reverse(files)[1:12]]


#=



ft = [Dates.Date(pstring(f)..., 1) for f in reverse(files)]
figure()
plot(ft, [length(r.SI) for r in recs], color = "green")
plot(ft, [length(r.SG) for r in recs], color = "blue")
plot(ft, [length(r.EG) for r in recs], color = "red")
ylim(0, 100)
=#

figure()
ax = subplot(;projection = ccrs.PlateCarree())
r = recs[10] # for March, 1980
lons = [];lats=[];mlds =[]
for k in keys(r)
    @show hasmld = findall(x->x==0, isnan.([r_.MLD for r_ in r[k]]))
    mldtemp = []
    for i in hasmld
        ind = findmin(abs.(r[k][i].MLD .- r[k][i].depth))[2]
        global s = scatter(r[k][i].lon, r[k][i].lat, c = r[k][i].T[ind], vmin = -2, vmax = 10)
    end
    
end
colorbar(s)
ax.coastlines()




#=

recnum = 1
figure()
for (i,k) in enumerate(keys(boxes))
    subplot(3, length(keys(boxes)), i)
    [plot(recs[recnum][k][i].T, recs[recnum][k][i].depth) for i in 1:length(recs[recnum][k])]
    title(k)
    ylim(1000,0)
    
    subplot(3, length(keys(boxes)), i + 3)
    [plot(recs[recnum][k][i].S, recs[recnum][k][i].depth) for i in 1:length(recs[recnum][k])]
    gca().invert_yaxis()
    ylim(1000,0)

    
    subplot(3, length(keys(boxes)), i + 6)
    for i in 1:length(recs[recnum][k])
        z = recs[recnum][k][i].depth
        nandices = (!).(isnan.(z))
        z = z[nandices]
        T= recs[recnum][k][i].T[nandices]
        S = recs[recnum][k][i].S[nandices]
        lat = recs[recnum][k][i].lat
        lon = recs[recnum][k][i].lon
        p = gsw_p_from_z.(-1*z, lat)
        SA = gsw_sa_from_sp.(S, p, lon, lat)
        CT = gsw_ct_from_t.(SA, T, p)
        ρ = gsw_rho.(SA, CT, p)
        try
            ml = hat.HolteAndTalley(p, T, S, ρ)
            println("success")
        catch
            println("nope")
        end
        
        plot(ρ,z)
    end
    
    gca().invert_yaxis()
    ylim(1000,0)
    
end
suptitle(ft[recnum])

figure()
subplot(;projection = ccrs.PlateCarree())
[scatter(r.lon, r.lat, color = "blue") for r in recs[recnum].SG]
[scatter(r.lon, r.lat, color = "green") for r in recs[recnum].SI]
[scatter(r.lon, r.lat, color = "red") for r in recs[recnum].EG]
gca().coastlines()
title(ft[recnum])

=#

#T and S diagrams for EN4 
#=
path = "/home/brynn/Code/oceanFP/data/EN4/analyses/"

wintermonths = collect(1:3)
function winterannual(x::Vector{T}, time) where T
    data = Vector{T}(undef, length(unique(year.(time))))
    for (i, y) in enumerate(unique(year.(time)))
        ind = intersect(findall(x->x == y, year.(time)), findall(x->x ∈ wintermonths, month.(time)))
        data[i] = mean(x[ind]) 
    end
    return data 
end

files = readdir(path)

nc = NCDataset(path * files[1]) 

lab_pt = (360-51, 57)
#lab_pt = (360-53, 61)

lon = nc["lon"][:]
lat = nc["lat"][:]
inds = [findall(x->x == lp, l)[1] for (lp, l) in zip(lab_pt, [lon, lat])]
close(nc) 
time = []
temp = []
salinity = []
for f in files
    nc = NCDataset(path * f) 
    push!(time, nc["time"][1])
    push!(temp, nc["temperature"][inds..., 1, 1])
    push!(salinity, nc["salinity"][inds..., 1,1])
    close(nc) 
end

using PyPlot
figure()
title("51°W, 57°N (Labrador Sea) EN4 JFM T-S") 
ind = findall(x->1920<x<1995, unique(year.(time)))
splot = winterannual(salinity, time)[ind]
tplot = winterannual(temp, time)[ind] .- 273.15
timeplot = unique(year.(time))[ind]
plot(splot, tplot, color = "black", alpha = 0.5)
s = scatter(splot, tplot, c = timeplot, cmap = "viridis", vmin = 1920, vmax = 1994)
[text(x,y,s) for (x,y,s) in zip(splot, tplot, string.(timeplot))]
colorbar(s)
xlabel("Mean Annual Salinity [g/kg]")
ylabel("Mean Annual Temperature [" * L"^\circ" * "C]")
savefig("/home/brynn/Code/OPTinv.jl/plots/en4tands.png") 

hadsst_path = "/home/brynn/Downloads/HadISST_sst.nc"
nc = NCDataset(hadsst_path)
lat = nc["latitude"][:]
lon = nc["longitude"][:]
lab_pt = (-51.5, 57.5)
inds = [findall(x->x == lp, l)[1] for (lp, l) in zip(lab_pt, [lon, lat])]
sst = nc["sst"][inds..., :]
wsst =  winterannual(sst, nc["time"][:])[begin:end-1]
timeplot_hadsst = unique(year.(nc["time"][:]))[begin:end-1]
close(nc)
figure()
plot(timeplot_hadsst, wsst, label = "HadISST")
plot(timeplot, tplot, label = "EN4")

title("51°W, 57°N (Labrador Sea) JFM Mean Temperature")

using CSV, DataFrames
df = CSV.read("/home/brynn/Code/OPTinv.jl/data/TLS1998Fig7Digitized.csv", DataFrame)
scatter(df[!, "Year"] .+ 1900, df[!, "Temperature"], label = "TLS1998", color = "green")
legend()
savefig("/home/brynn/Code/OPTinv.jl/plots/labtemp.png") 
=#


#Code Snippet to Get 1 MLD 
#=
using NaNMath
depths = [NaNMath.maximum(convert(Vector{Float32}, recs[1].SG[i].depth)) for i in 1:length(recs[1].SG)]
findall(x->x>1000, depths)
i = 69;figure();plot(recs[1].SG[i].T, recs[1].SG[i].depth);gca().invert_yaxis();ylabel("Depth [m]");xlabel("Temperature")

p = gsw_p_from_z.(-1 * recs[1].SG[i].depth, recs[1].SG[i].lat)
T = recs[1].SG[i].T
h = hat.HolteAndTalley(p[(!).(isnan.(p))], T[(!).(isnan.(T))])
h.tempMLD
=#
=#
