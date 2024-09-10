#=
Generate a Hovmoller plot using EN4 reanalyses for a point in the Labrador Sea
=#
import Pkg;Pkg.activate("../../") 
using NCDatasets, PythonPlot, OPTinv, Dates, Statistics, DrWatson, CSV, DataFrames

path = "/home/brynn/Code/oceanFP/data/EN4/analyses"
latp = 60.164480
lonp = 360-54.655807
nct = NCDataset(joinpath(path, readdir(path)[1]))
lat = nct["lat"][:] #degN
lon = nct["lon"][:] #degE
latpi = findmin(abs.(lat .- latp))[2]
lonpi = findmin(abs.(lon .- lonp))[2]
templatet = nct["temperature"][:, :, 1, 1]
templatet[ismissing.(templatet)] .= NaN
templatet = convert(Matrix{Float64}, templatet')
figure();contourf(lon, lat, templatet)
scatter(lon[lonpi], lat[latpi], color = "red")
dataind = (!).(ismissing.(nct["temperature"][lonpi, latpi, :, 1]))
depth = nct["depth"][dataind] #m 
nct["temperature"][lonpi, latpi, dataind, 1]
mat = Array{Float32}(undef, length(depth), length(readdir(path)), 2)
time = Vector{DateTime}(undef, length(readdir(path)))
close(nct)
for (i,f) in enumerate(readdir(path))
    println("Opening file " * string(i) * " out of " * string(length(readdir(path))))
    nc = NCDataset(joinpath(path, f)) 
    time[i] = nc["time"][1]
    mat[:, i, 1] = nc["temperature"][lonpi, latpi, dataind, 1]
    mat[:, i, 2] = nc["salinity"][lonpi, latpi, dataind, 1]
    close(nc)
end
size(mat)
t = unique(year.(time))[begin:end-1]
matann = hcat([mean(mat[:, (i-1)*12+1:i*12, :], dims = 2)[:, :,:] for i in 1:length(t)]...)

tdec = 1900:10:2010
matdec = hcat([mean(mat[:, (i-1)*12*10+1:i*12*10, :], dims = 2)[:, :,:] for i in 1:length(tdec)]...)

figure(figsize = (8,8))
ax = subplot(211)
cf = contourf(t, depth, matann[:, :, 1] .- mean(matann[:, :, 1], dims = 2), cmap = "coolwarm", vmin = -3, vmax = 3)
  cb = colorbar(cf)
cb.set_label("Temperature Anomaly [°C]") 
c = contour(t, depth, matann[:, :, 1] .- mean(matann[:, :, 1], dims = 2), colors = "black", vmin = -3, vmax = 3, linewidths = 0.5)
#ax.clabel(c)
ylim(3250,0)
ax = subplot(212)
cf = contourf(t, depth, matann[:, :, 2] .- mean(matann[:, :, 2], dims = 2), cmap = "coolwarm", vmin = -0.6, vmax = 0.6)
c = contour(t, depth, matann[:, :, 2] .- mean(matann[:, :, 2], dims = 2), colors = "black", linewidths = 0.5, vmin = -0.6, vmax = 0.6)
cb = colorbar(cf)
cb.set_label("Salinity Anomaly [psu]")
ylim(3250,0)
tight_layout()

figure(figsize = (8,8))
lev = 25
ax = subplot(211)
cf = contourf(tdec, depth, matdec[:, :, 1] .- mean(matdec[:, :, 1], dims = 2), cmap = "coolwarm", vmin = -1.5, vmax = 1.5, levels = lev)
cb = colorbar(cf)
cb.set_label("Temperature Anomaly [°C]") 
c = contour(tdec, depth, matdec[:, :, 1] .- mean(matdec[:, :, 1], dims = 2), colors = "black", vmin = -1.5, vmax = 1.5, linewidths = 0.5,levels = lev)
ylabel("Depth [m]") 
#ax.clabel(c)
ylim(3250,0)
ax = subplot(212)
cf = contourf(tdec, depth, matdec[:, :, 2] .- mean(matdec[:, :, 2], dims = 2), cmap = "coolwarm", vmin = -0.32, vmax = 0.32,levels = lev)
c = contour(tdec, depth, matdec[:, :, 2] .- mean(matdec[:, :, 2], dims = 2), colors = "black", linewidths = 0.5, vmin = -0.32, vmax = 0.32,levels = lev)
cb = colorbar(cf)
cb.set_label("Salinity Anomaly [psu]")
xlabel("Time [years CE]")
ylabel("Depth [m]") 

#ax.clabel(c)
ylim(3250,0)
tight_layout()
savefig(plotsdir("labrador_hovmoller.png"))

figure()
s = matann[findmin(abs.(depth .- 250))[2], :, :]
x = s[:, 2] # .- mean(s[:, 2])
y = s[:, 1] .- 273.15# .- mean(s[:, 1])
scat = scatter(x,y, c = t, cmap = "Reds", s = 100)
#plot(x,y, color = "gray", zorder = 0)
#[text(x = xi, y = yi, s = ti, fontsize = 12) for (xi, yi, ti) in zip(x,y,tdec)]
cb = colorbar(scat)
cb.set_label("Time [year CE]", fontsize = 15)
xlabel("Salinity Anomaly [psu]", fontsize = 15)
ylabel("Temperature Anomaly [°C]", fontsize = 15)
title("EN4: 60.16°N, 54.66°W, 250m", fontsize = 20)
savefig(plotsdir("en4_t_and_s_250m.png"))

figure()
title("EN4: 60.16°N, 54.66°W, 1500m", fontsize = 20)
s = matann[findmin(abs.(depth .- 1500))[2], :, :]
t_tls = 1900 .+ df[!, "Year"]
tind = findall(x -> minimum(t_tls) < x < maximum(t_tls), t)
x = s[tind, 2] # .- mean(s[:, 2])
y = s[tind, 1] .- 273.15# .- mean(s[:, 1])
scat = scatter(x,y, c = t[tind], cmap = "Reds", s = 100, vmin = minimum(t_tls), vmax = maximum(t_tls), edgecolors = "black", label = "EN4")
cb = colorbar(scat)
cb.set_label("Time [year CE]", fontsize = 15) 
xlabel("Salinity Anomaly [psu]", fontsize = 15)
ylabel("Temperature Anomaly [°C]", fontsize = 15) 
df = CSV.read(OPTinv.datadir("TLS1998Fig7Digitized.csv"), DataFrame)
scatter(df[!, "Salinity"], df[!, "Temperature"], c = t_tls, cmap = "Reds", marker = "*", edgecolors = "black", s = 150, vmin = minimum(t_tls), vmax = maximum(t_tls), label = "TLS1998")
legend()
savefig(plotsdir("en4_t_and_s_1500m.png"))

