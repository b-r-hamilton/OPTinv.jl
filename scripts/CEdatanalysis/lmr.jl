import Pkg
Pkg.activate("../../")
using NCDatasets, PythonPlot, DrWatson, Statistics, PythonCall, DateFormats, RollingFunctions, Dates, OPTinv, PaleoData

#https://www.ncei.noaa.gov/pub/data/paleo/reconstructions/tardif2019lmr/v2_0/sst_MCruns_ensemble_mean_LMRv2.0.nc
lmrdataset = loadLMR("sst")
time = lmrdataset["time"][:]; lon = lmrdataset["lon"][:]; lat = lmrdataset["lat"][:]
if ! @isdefined sst
    sst = mean(makeNaN(lmrdataset["sst"][:, :, :, :]), dims = 3)[:, :, 1, :]
end

years = 1500:50:1900
figure(figsize = (10, 8))
for (i, y) in enumerate(years)
    ax = subplot(3,3,i, projection = ccrs.PlateCarree())
    title(string(y) *"-" * string(y+50))
    ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "110m", edgecolor="k", facecolor="gray"))
    ax.set_extent([-80,30, 0, 90])
    plotme = mean(sst[:, :, y+50], dims = 3)[:, :, 1]'
    plotme, plotlon = cu.add_cyclic_point(plotme, coord = lon)
    cf = contourf(plotlon, lat, plotme, cmap = cm.balance, levels = -1:0.1:1)
    c = contour(plotlon, lat, plotme, colors = "black", levels = -1:0.1:1)#, vmin = -vm, vmax = vm, levels = lev )
    ax.clabel(c) 
    #contour(lon
    #colorbar(cf)
end
tight_layout()
savefig(plotsdir("lmrmaps.png"))
μ = [mean(sst[:, :, i][(!).(isnan.(sst[:, :, i]))]) for i in 1:length(time)]

μbox = boxmean(sst, lat, lon, 50, 90, 360-50, 20) 


#https://www.wdc-climate.de/ui/q?hierarchy_steps_ss=ModE-RA_s14203-18501&entry_type_s=Dataset
nc = NCDataset(OPTinv.datadir("ModE-RA_ensmean_temp2_anom_wrt_1901-2000_1421-2008_mon.nc"))

temp = nc["temp2"][:, :, :]
time = nc["time"][:]
ind1980 = findall(x->Dates.DateTime(1980,1,1) < x < Dates.DateTime(1981,1,1), time)
temp .-= mean(temp[:, :, ind1980], dims =3)
lon = nc["longitude"][:]
lat = nc["latitude"][:]

figure(figsize = (10, 8))
for (i, y) in enumerate(years) 
    time_ind = findall(x->x == 1, Dates.DateTime(y,1,1) .< time .< Dates.DateTime(y+50, 1, 1))
    ax = subplot(3,3,i, projection = ccrs.PlateCarree())
    ax.set_extent([-80,30, 0, 90])

    title(string(y) *"-" * string(y+50))
    #ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "110m", edgecolor="gray"))
    ax.coastlines(color = "gray")
    plotme = mean(temp[:, :, time_ind], dims = 3)[:, :, 1]'
    plotme, plotlon = cu.add_cyclic_point(plotme, coord = lon)
    cf = contourf(plotlon, lat, plotme, cmap = cm.balance, levels = -1:0.2:1)
    c = contour(plotlon, lat, plotme, colors = "black", levels = -1:0.2:1)
    ax.clabel(c) 
end
tight_layout()
savefig(plotsdir("moderamaps.png"))
roll = 50
figure()
μmodera = [mean(temp[:, :, i][(!).(isnan.(temp[:, :, i]))]) for i in 1:length(time)]
ydtime = yeardecimal.(time)
ydtimer = rolling(mean, ydtime, roll * 12)
μmoderar = rolling(mean, μmodera, roll * 12)
plot(yeardecimal.(time), μmodera, color = "blue", alpha = 0.5)
plot(ydtimer, μmoderar, color = "blue",label ="ModE-RA")



plot(collect(-1:1:1999), μ, color = "black", alpha = 0.5)
length(μ)
length(-1:1:2000)

plot(collect(-1:1:1999), μbox, color = "red", alpha = 0.5)


μr = rolling(mean, μ, roll)
μboxr = rolling(mean, μbox, roll)
timer = rolling(mean, collect(-1:1:1999), roll)
plot(timer, μr, color = "black",label = "LMR: Global")
plot(timer, μboxr, color = "red", label = "LMR: 50°W-20°E, 50°N-90°N")
xlabel("Time [years CE]")
ylabel("Mean Temperature Anomaly From 1900-2000 Mean [K]")
xlim(1000,2000)
legend()
savefig(plotsdir("reanalysis_timeseries.png"))
