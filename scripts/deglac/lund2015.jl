import Pkg
Pkg.activate("../../")
using OPTinv, PaleoData, LinearAlgebra, TMI, PythonPlot, Interpolations
locs, data = loadLund2015()

#TMIversion = "modern_180x90x33_GH11_GH12"
#TMIversion = versionlist()[1]
versionlist()
TMIversions = ["modern_90x45x33_GH10_GH12", "LGM_90x45x33_G14", "LGM_90x45x33_G14A", "LGM_90x45x33_GPLS1", "LGM_90x45x33_GPLS2", "LGM_90x45x33_OG18"]
struct TMIdat
    smat::SVD
    pv::Vector{Float64} 
    γ::Grid
    obsSH::Vector{Float64}
    obsSSH::Vector{Float64}
    obsNH::Vector{Float64}
end
TMIdats = Vector{TMIdat}(undef, length(TMIversions))
for (i, TMIversion) in enumerate(TMIversions)
    println(TMIversion)
    mat, smat = OPTinv.generatemodes(locs, TMIversion = TMIversion, func = svd)
    pv = smat.S.^2 ./ sum(smat.S.^2)
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

    sh = γ.wet[:,:,1]
    nh = γ.wet[:,:,1]
    ssh = γ.wet[:,:,1]
    sh[:, γ.lat .> 0] .= 0
    nh[:, γ.lat .< 0] .= 0
    ssh[:, γ.lat .> -60] .= 0
    SH = BoundaryCondition(sh, γ.axes[1:2], 0.0, 3,1, γ.wet[:, :,1], :SH, "southern hemi", "southern hemi")
    NH = BoundaryCondition(nh, γ.axes[1:2], 0.0, 3,1, γ.wet[:, :,1], :NH, "northern hemi", "northern hemi")
    sSH = BoundaryCondition(ssh, γ.axes[1:2], 0.0, 3,1, γ.wet[:, :,1], :SSH, "southern southern hemi", "southern southern hemi")
    siSH = steadyinversion(Alu, SH, γ)
    siSSH = steadyinversion(Alu, sSH, γ)
    siNH = steadyinversion(Alu, NH, γ)
    N = length(locs)
    wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,N)
    [wis[i] = interpindex(locs[i],γ) for i in 1:N]
    obsSH = observe(siSH, wis, γ)
    obsSSH = observe(siSSH, wis, γ)
    obsNH = observe(siNH, wis, γ)
    TMIdats[i] = TMIdat(smat, pv, γ, obsSH, obsSSH, obsNH) 
end

#=
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
    ylim(350, 4100)
    xlim(-1,1)                      
    ax.invert_yaxis()
end
=#
depths = [l[3] for l in locs]
for (TMIv, TMId) in zip(TMIversions, TMIdats)
    println(TMIv) 
    println("First three modes capture " * string(round(sum(TMId.pv[begin:3]), sigdigits = 4) * 100) * "% of variance")

    proj = ccrs.PlateCarree()
    figure(figsize = (12,4))
    suptitle(TMIv) 
    for i in 1:3
        ax = subplot(1,3,i, projection = proj)
        ax.coastlines()
        title("Mode " * string(i) * ": " * string(round(TMId.pv[i] * 100, sigdigits = 3)) * "% var. exp.")
        plotme = γreshape(TMId.smat.Vt[i, :], TMId.γ)'
        plotme, lon_ = cu.add_cyclic_point(plotme, coord = TMId.γ.lon)
        cf = contourf(lon_, TMId.γ.lat, plotme, transform = proj)
        colorbar(cf, fraction = 0.02)
    end
    tight_layout()
    savefig("spatialmodes" * TMIv *".png")
    
    figure(figsize = (4,4))
    suptitle(TMIv) 
    for i in 1:3
        ax = subplot(1,3, i)
        scatter(TMId.smat.U[:, i], depths, c = TMId.smat.U[:, i], cmap = cm.balance, vmin = -1, vmax = 1, edgecolors = "gray", 100)
        yt = ax.get_yticks()
        if i != 1
            ax.set_yticks(yt, fill("", length(yt)))
        end
        title("Mode " * string(i))
        ylim(350, 4100)
        xlim(-1,1)                      
        ax.invert_yaxis()
    end
    savefig("pv" * TMIv * ".png")
end

linestyles = ["solid", (0, (1,1)), "dotted", "dashed", "dashdot", (0,(1,10))]
figure(figsize = (8,8))
for (linestyle, TMIv, TMId) in zip(linestyles, TMIversions, TMIdats)
    println(TMIv)
    plot(TMId.obsNH[sortperm(depths)], sort(depths), label = "NADW: " * TMIv, color = "tab:blue", linestyle = linestyle)
    plot(TMId.obsSSH[sortperm(depths)], sort(depths), label = "AABW: " * TMIv, color = "tab:green", linestyle = linestyle)
    plot((TMId.obsSH .- TMId.obsSSH)[sortperm(depths)], sort(depths), label = "AAIW: " * TMIv, color = "tab:orange", linestyle = linestyle)
end

errorbar([0.7,0.56,0.42], [1820,2082,2296], xerr = [0.22, 0.20, 0.19], label = "LGM (from d18O, Lund2015)", color = "darkblue", fmt = "x", capsize = 5)
errorbar([0.52,0.45,0.35], [1820,2082,2296], xerr = [0.1,0.1,0.09], label = "HS1 (from d18O, Lund2015)", color = "cornflowerblue", fmt = "^", capsize = 5)
ylabel("Depth [m]")
#legend()
gca().invert_yaxis()
xlabel("Proportion Water")
legend(loc="center left", bbox_to_anchor=(1, 0.5))
tight_layout()
savefig("TMIvLundpercent.png")

re2005 = loadRickabyandElderfield2005()
function cleanup(v)
    v[typeof.(v) .!= Float64] .= NaN
    return convert(Vector{Float64}, v)
end
reage = cleanup(re2005.age)
red18O = cleanup(re2005.d18o)
red13C = cleanup(re2005.d13c)
using NaNMath, Statistics
NADWholo = NaNMath.mean(red18O[findall(x->0<x<7, reage)])
AAIWholo = mean(data.ggc14.d18Ocwuell250[findall(x->0<x<7, data.ggc14.age_calkaBP)])
AABWholo =mean(data.ggc22.d18Ocwuell250[findall(x->0<x<7, data.ggc22.age_calkaBP)])

figure()
holod18O = [mean(data[k].d18Ocwuell250[findall(x->0<x<7, data[k].age_calkaBP)]) for k in keys(data)]
scatter(holod18O, depths, color = "black", label = "Mean Holocene")
#scatter(lgmd18O, depths, color = "black", marker = "x")
gca().invert_yaxis()
for (linestyle, TMIv, TMId) in zip(linestyles, TMIversions, TMIdats)
    holod18Orec = @. TMId.obsSSH * AABWholo + (TMId.obsSH .- TMId.obsSSH) * AAIWholo + TMId.obsNH * NADWholo
    holod18Orec2end = @. TMId.obsSH * AABWholo + TMId.obsNH * NADWholo
    plot(holod18Orec[sortperm(depths)], sort(depths), color = "red", linestyle = linestyle, label = "3: " * TMIv)
    plot(holod18Orec2end[sortperm(depths)], sort(depths), color = "pink", linestyle = linestyle, label = "2: " * TMIv)
end
legend()
ylabel("Depth [m]")
xlabel("d18Oc [‰ VPDB]")
savefig("holorec.png") 

lgmd18O = [mean(data[k].d18Ocwuell250[findall(x->19<x<22, data[k].age_calkaBP)]) for k in keys(data)]
lgmd13C = [mean(data[k].d13Ccwuell250[findall(x->19<x<22, data[k].age_calkaBP)]) for k in keys(data)]
NADWlgmd18O = NaNMath.mean(red18O[findall(x->19<x<22, reage)])
AAIWlgmd18O = mean(data.ggc14.d18Ocwuell250[findall(x->19<x<22, data.ggc14.age_calkaBP)])
AABWlgmd18O =mean(data.ggc22.d18Ocwuell250[findall(x->19<x<22, data.ggc22.age_calkaBP)])

NADWlgmd13C = NaNMath.mean(red13C[findall(x->19<x<22, reage)])
AAIWlgmd13C = mean(data.ggc14.d13Ccwuell250[findall(x->19<x<22, data.ggc14.age_calkaBP)])
AABWlgmd13C =mean(data.ggc22.d13Ccwuell250[findall(x->19<x<22, data.ggc22.age_calkaBP)])

mat = [NADWlgmd18O AAIWlgmd18O AABWlgmd18O
       NADWlgmd13C AAIWlgmd13C AABWlgmd13C
       1 1 1]

using NonNegLeastSquares
res = hcat([nonneg_lsq(mat, [lgmd18O[i], lgmd13C[i], 1]) for i in 1:12]...)
#how well did we fit the obs? 
[mat * res[:, i] for i in 1:12] .- [[lgmd18O[i], lgmd13C[i], 1] for i in 1:12]

figure()
subplot(1,3,1)
plot(res[1,:][sortperm(depths)], sort(depths), label = "%NADW", color = "tab:blue")
for (linestyle, TMIv, TMId) in zip(linestyles, TMIversions, TMIdats)
    plot(TMId.obsNH[sortperm(depths)], sort(depths), alpha = 0.5, color = "tab:blue", label = "NADW: " * TMIv, linestyle = linestyle)
end
ylabel("depth [m]")
xlabel("proportion of contribution")
title("NADW")
errorbar([0.7,0.56,0.42], [1820,2082,2296], xerr = [0.22, 0.20, 0.19], label = "LGM (from d18O, Lund2015)", color = "darkblue", fmt = "x", capsize = 5)
xlim(-0.01,1)
gca().invert_yaxis()
subplot(1,3,2)

plot(res[2,:][sortperm(depths)], sort(depths), label = "%AAIW", color = "tab:orange")
for (linestyle, TMIv, TMId) in zip(linestyles, TMIversions, TMIdats)
    plot((TMId.obsSH .- TMId.obsSSH)[sortperm(depths)], sort(depths), alpha = 0.5, color = "tab:orange", label = "AAIW: " * TMIv, linestyle = linestyle)
end
xlabel("proportion of contribution")
title("AAIW")
xlim(-0.01,1)
gca().invert_yaxis()
subplot(1,3,3)
for (linestyle, TMIv, TMId) in zip(linestyles, TMIversions, TMIdats)
    plot(TMId.obsSSH[sortperm(depths)], sort(depths), alpha = 0.5, color = "tab:green", label = "AABW: " * TMIv, linestyle = linestyle)
end
plot(res[3,:][sortperm(depths)], sort(depths), label = "%AABW", color = "tab:green")

xlabel("proportion of contribution")
title("AABW")
xlim(-0.01,1)
gca().invert_yaxis()
#legend()
tight_layout()
savefig("lgmprop.png")
