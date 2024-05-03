import Pkg
Pkg.activate("../")


using OPTinv, Unitful, DimensionalData, BLUEs, Statistics, UnitfulLinearAlgebra, LinearAlgebra, PyPlot, Revise, DrWatson, PyCall, TMI, CSV, DataFrames, JLD2, NaNMath, Dates, RollingFunctions
import Measurements.value as value
import OPTinv.Est


ccrs = pyimport("cartopy.crs")
cm = pyimport("cmocean.cm")
cfeature = pyimport("cartopy.feature")
cu = pyimport("cartopy.util")

core_list = Symbol.("MC".* string.([28, 26, 25, 22, 21, 20, 19, 10, 9,13, 14]) .* "A")
#load in ℳ and surface spatial modes 
filename = "nmf.jld2"
ℳnmf, spatialmodesnmf = loadM(filename, core_list)
filename = "svd.jld2"
ℳsvd, spatialmodessvd = loadM(filename, core_list)

if ! @isdefined A 
    TMIversion = "modern_180x90x33_GH11_GH12"
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
end


lon = ustrip.(γ.lon); lat = ustrip.(γ.lat)
#lev = -0.05:0.01:0.05
for (sm, t) in zip([spatialmodessvd, spatialmodesnmf], ["SVD", "NMF"]) 
    figure(figsize = (6, 6))
    for i in 1:11 
        ax = subplot(3,4,i, projection = ccrs.PlateCarree())
        plotme = γreshape(sm[i, :], γ)'
        plotme, lon = cu.add_cyclic_point(plotme, coord = ustrip.(γ.lon))
        pm = plotme[(!).(isnan.(plotme))]
        vm = maximum(abs.([minimum(pm), maximum(pm)]))
        cf = contourf(lon, lat, plotme, cmap = cm.balance, vmin = -vm, vmax = vm)
        #colorbar(cf)
        
        #c = contour(lon, lat, plotme, colors = "black", levels = lev)
        #ax.clabel(c)
        ax.set_extent([-80, 30, -90, 90])
        ax.coastlines(color = "gray")
        title(i)
    end
    tight_layout()
    savefig(plotsdir(t * "_modes.png"))
end


        
                 
