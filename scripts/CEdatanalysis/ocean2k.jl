#=
Plots every Ocean2k record
Computes dSST/dt for each record in ROI and plots on map 
=#
import Pkg;Pkg.activate("../")
using OPTinv, Statistics, DataFrames, PythonPlot, PaleoData, NaNMath, DrWatson 

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


figure()
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
    post1000 = findall(x->x>1000, age)
    if NaNMath.maximum(age) > 1500
        lls, C = linearleastsquares(age[post1000], sst[post1000])
        @show lls[1] * 1000
        d[r] = (lat_, lon_, lls[1] * 1000)
    end
    
    plot(age, sst, label = lab)
end
xlabel("Time [years CE]")
ylabel("SST [°C]")
legend()
savefig(plotsdir("ocean2k_recs.png"))
                   
figure()
subplot(projection = ccrs.PlateCarree())
for (i, k) in enumerate(keys(d))
    s = scatter(d[k][2], d[k][1], c = d[k][3], vmin = -3, vmax = 3, cmap = cm.balance, s = 100)
    if i == 1
        cb = colorbar(s)
        cb.set_label("dSST/dt [°C/kyr]")
    end
    
end
gca().coastlines()
savefig(plotsdir("ocean2k.png"))
