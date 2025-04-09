#=
Script to generate σ vector of mode magnitude uncertainties from CESM output
Also makes supplementary plot
=# 
#cd("/home/brynn/Code/OPTinv.jl/scripts/CEdatanalysis")
using NCDatasets, OPTinv, Interpolations, TMI, PythonPlot, DrWatson, JLD2, Statistics, PythonCall
"""
function ginterp

    interpolate z(x,z) to (x2,y2) grid
"""
function ginterp(x,y,z::Array{Float32, 2}, x2, y2)
    x2 = matchvec(x2, x)
    y2 = matchvec(y2, y) 
    itp = LinearInterpolation((x, y), z)
    z2 = [itp(x,y) for x in x2, y in y2]
    return z2, x2, y2
end
"""
matchvec: trim a vector x1 so that it does not exceed extrema of x2
"""
matchvec(x1, x2) =  x1[x1 .< maximum(x2) .&& x1 .> minimum(x2)]

#read CESM2 SST data from realization `r11i1p1f1`
urls = ["https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/NCAR/CESM2/historical/r11i1p1f1/Omon/thetao/gr/v20190514/thetao_Omon_CESM2_historical_r11i1p1f1_gr_185001-189912.nc", 
        "https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/NCAR/CESM2/historical/r11i1p1f1/Omon/thetao/gr/v20190514/thetao_Omon_CESM2_historical_r11i1p1f1_gr_190001-194912.nc",
        "https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/NCAR/CESM2/historical/r11i1p1f1/Omon/thetao/gr/v20190514/thetao_Omon_CESM2_historical_r11i1p1f1_gr_195001-199912.nc",
        "https://esgf.ceda.ac.uk/thredds/dodsC/esg_cmip6/CMIP6/CMIP/NCAR/CESM2/historical/r11i1p1f1/Omon/thetao/gr/v20190514/thetao_Omon_CESM2_historical_r11i1p1f1_gr_200001-201412.nc"]        
@time thetas = [NCDataset(url)["thetao"][:, :, 1, :] for url in urls]
thetas = cat(thetas..., dims = 3)
thetas[ismissing.(thetas)] .= NaN
thetas = convert(Array{Float32, 3}, thetas)

lat = NCDataset(urls[1])["lat"][:]
lon = NCDataset(urls[1])["lon"][:]
times = vcat([NCDataset(url)["time"][:] for url in urls]...)
interval = 120
thetadec = cat([mean(thetas[:, :, ((i-1)*interval)+1:i*interval], dims = 3) for i in 1:convert(Int64, floor(length(times) / interval))]..., dims = 3)
timesdec = collect(1850:interval/12:2000)

corenums_full = Symbol.("MC" .* string.([28, 26, 25, 22, 21, 20, 19, 10, 9,13,14]) .* "A")
filename = "svd.jld2"
cd("../../src")
ℳ, spatialmodes = loadM(filename, corenums_full) #V
N = size(spatialmodes)[2]
spatialmodesglobal = reshape(fill(sqrt(1/N), N), (1, N)) #Vglobal 

TMIversion = "modern_180x90x33_GH11_GH12"
γ = Grid(download_ncfile(TMIversion))
tmodes = Matrix{Float64}(undef, size(spatialmodes)[1], length(timesdec))
tmodesglobal = Matrix{Float64}(undef, size(spatialmodesglobal)[1], length(timesdec))

_, newlon, newlat = ginterp(lon, lat, thetadec[:, :, 1], γ.lon, γ.lat)
for i in 1:length(timesdec)
    mat, newlon, newlat = ginterp(lon, lat, thetadec[:, :, i], γ.lon, γ.lat)
    nind = (!).(isnan.(mat[γ.wet[:, :, 1]]))
    tmodes[:, i] .= spatialmodes[:, nind]' \ mat[γ.wet[:, :, 1]][nind]
    tmodesglobal[:, i] .= spatialmodesglobal[:, nind]' \ mat[γ.wet[:, :, 1]][nind]
end
svdmags = tmodes'
jldsave(DrWatson.datadir("modemags.jld2"); svdmags)
svdmags = tmodesglobal'
jldsave(DrWatson.datadir("modemags_global.jld2"); svdmags)

jld = jldopen(DrWatson.datadir("modemags_div5.jld2"))
mags = jld["svdmags"]    
σθ = vec(std(mags, dims = 1)) * K
println("dividing σθ by 5") 
σθ ./= 5

figure(figsize = (10,4))
ax = subplot(1,2,1)
ax.plot(timesdec, tmodes[1, :], color = "black", linewidth = 3)
ax.set_ylabel(L"u_\mathrm{model}(t) []", fontsize = 15)
yt = ax.get_yticks()
yticks(yt, fontsize = 12)
ax1 = twinx()
#ax1.plot(timesdec, [boxmean(γreshape(spatialmodes[1,:] .* tmodes[1,i], γ), newlat, newlon, 60,80, 0, 20) for i in 1:length(timesdec)], color = "red", linewidth = 2, linestyle = "dashed")
ax1.plot(timesdec, boxmean(thetadec, lat, lon, 60, 80, 0, 20), linewidth = 2, color = "red")
ax.invert_yaxis()
yt = 4.4:0.2:6
ax1.set_yticks(yt, yt, fontsize = 12, color = "red")
ax1.set_ylabel("T [°C]", fontsize = 15, color = "red")
xt = pyconvert(Vector{Int64}, ax1.get_xticks())
ax.set_xticks(xt,xt, fontsize = 12)
ax.set_xlabel("Time [years CE]", fontsize = 15)
text(x = 2005, y = 5.8, s = "A", fontweight = "bold", fontsize = 30)

ax3 = subplot(1,2,2);plot(1:11, std(tmodes, dims = 2), color = "black", ".-", markersize = 12)
ax3.set_xticks(1:11, 1:11, fontsize = 12)
ax3.set_xlabel("Mode Number", fontsize = 15)
ax3.grid()
ax3.set_yticks(1:6, 1:6, fontsize = 12)
ax3.set_ylabel("Standard Deviation of "* L"u_\mathrm{model}" *" []", fontsize = 15)
ax3.text(x = 10.5, y = 5.95, s = "B", fontweight = "bold", fontsize = 30)
tight_layout()
savefig(DrWatson.plotsdir("sigma.png"), dpi = 600)
