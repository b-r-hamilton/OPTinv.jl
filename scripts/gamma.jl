#=
Script to explore sensitivity of theoretical ΔT equation to γ
=# 
import Pkg;Pkg.activate("../")
using PythonPlot
Δδ¹⁸O = 0.08
α = -0.224
ΔT(γ) = Δδ¹⁸O/(α + 1/γ)
γcrit = -1/α
γ1 = range(start = 0, stop = γcrit-0.01, length = 100)
γ2 = range(start = γcrit+0.01, stop = 40, length = 100)
figure(figsize = (8,4))
subplot(1,2,1) 
plot(γ1, ΔT.(γ1), color = "black", linewidth = 3)
plot(γ2, ΔT.(γ2), color = "black", linewidth = 3)
vlines(x = γcrit, ymin = -2, ymax = 2, color = "gray", linestyle = "dotted", linewidth = 2)
ylim(-2,2)
title(L"\Delta T = \frac{\Delta \delta^{18}\mathrm{O}_\mathrm{seawater}}{\alpha + \gamma^{-1}}", fontsize = 20)
ylabel("ΔT [K]", fontsize = 15)
xlabel("γ [K ‰⁻¹]", fontsize = 15)
yt = range(start = -2, stop = 2, step = 0.5)
yticks(yt, yt, fontsize = 12)
xt = range(start = 0, stop = 40, step = 5)
hlines(y = Δδ¹⁸O/α, xmin = 0, xmax = 40, linestyle = "dashed", color = "gray", linewidth = 2)
xticks(xt, xt, fontsize = 12)
xlim(0,40)
poi = [14.93,35.01]
scatter(x = poi, y = ΔT.(poi), marker = "o", color = "white", s = 100, edgecolor = "black",zorder = 10000)
text(x = 35, y = 1.4, s = "A", fontsize = 30, weight = "bold")
subplot(1,2,2) 
γT = range(5, 20, length = 200)
γs = range(2, 3, length = 200)
mat = zeros(length(γT), length(γs))
for (i, γT_) in enumerate(γT)
    for (j, γs_) in enumerate(γs)
        mat[i,j] = ΔT(γT_ * γs_)
    end
end

c = contour(γT, γs, mat, cmap = "seismic",vmin = -0.7, vmax = 0.7, linewidths = 3)
clabel(c, fontsize = 12)
γspoi = [2.796, 2.257]
γTpoi = [5.34, 15.51]
scatter(γTpoi, γspoi, marker = "o", color = "white", s = 100, edgecolor = "black",zorder = 10000)
ylabel("γₛ [g kg⁻¹ ‰]", fontsize = 15)
xlabel("γₜ [K (g kg⁻¹)⁻¹]", fontsize = 15)
title(L"\Delta T = \frac{\Delta \delta^{18}\mathrm{O}_\mathrm{seawater}}{\alpha + (\gamma_s \gamma_t)^{-1}}", fontsize = 20)
yt = range(start = 2, stop = 3, step = 0.2)
yticks(yt, yt, fontsize = 12)
xt = range(start = 5, stop = 20, step = 5)
xticks(xt, xt, fontsize = 12)
text(x = 18, y = 2.85, s = "B", fontsize = 30, weight = "bold")
tight_layout()
savefig("../plots/gamma.png", dpi = 600)
