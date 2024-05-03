using NCDatasets, Dates, CSV, DataFrames, Statistics
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
