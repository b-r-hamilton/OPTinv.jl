import Pkg
Pkg.activate("../../")
using OPTinv, PaleoData, LinearAlgebra, TMI, PythonPlot

locs = loadLund2015()

mat, smat = OPTinv.generatemodes(locs, func = svd)

TMIversion = "modern_180x90x33_GH11_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

pv = smat.S.^2 ./ sum(smat.S.^2)
figure(figsize = (15,15))
for i in 1:12 
    subplot(4,3,i, projection = ccrs.PlateCarree())
    title(string(round(pv[i] * 100, sigdigits = 3)) * "% var. exp.")
    cf = contourf(γ.lon, γ.lat, γreshape(smat.Vt[i, :], γ)')
    colorbar(cf) 
end

depths = [l[3] for l in locs]
figure(figsize = (10,4))
for i in 1:12
    ax = subplot(1,12, i)
    scatter(smat.U[:, i], depths, c = smat.U[:, i], cmap = cm.balance, vmin = -1, vmax = 1, edgecolors = "gray", 100)
    yt = ax.get_yticks()
    if i != 1
        ax.set_yticks(yt, fill("", length(yt)))
    end
    xlim(-1,1)                      
    ax.invert_yaxis()
end


