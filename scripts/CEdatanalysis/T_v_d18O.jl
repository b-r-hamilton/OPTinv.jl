
import Pkg;Pkg.activate("../../")
using TMI, NCDatasets, PythonPlot, OPTinv, UnitfulLinearAlgebra, PythonCall, Statistics, DrWatson

TMIversion = "modern_180x90x33_GH11_GH12"
if ! @isdefined A
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
end


nc = NCDataset(TMIfile) 

d18O = nc["δ¹⁸Ow"][:, :, 1]
σd18O = nc["σδ¹⁸Ow"][:, :, 1]
T = nc["θ"][:, :, 1]
S = nc["Sp"][:, :, 1]

pairs = [[S, T], [d18O, S], [d18O, T]]
xlabels = ["S [psu]", "δ¹⁸O [‰]", "δ¹⁸O [‰]"]
ylabels = ["T [°C]", "S [psu]", "T [°C]"]
vy = ["T", "S", "T"]
vx = ["S", "δ¹⁸O", "δ¹⁸O"]
letters = ["A", "B", "C"]

figure(figsize = (13, 4))
for (i, (p, xlab, ylab, varx, vary, letter)) in enumerate(zip(pairs, xlabels, ylabels, vx, vy, letters))
    ax = subplot(1,3,i) 
    scatter(p[1], p[2], color = "black", alpha = 0.1)
    dind = (!).(isnan.(vec(p[1])))
    x = vec(p[1])[dind]
    y = vec(p[2])[dind]
    lls, llsc = linearleastsquares(x,y)
    ỹ = x.*lls[1] .+ lls[2]
    redind = findall(x->x>minimum(y), ỹ)

    plot(x[redind], ỹ[redind], color = "red", linewidth = 5)
    ylim(minimum(y), maximum(y))
    xlim(minimum(x), maximum(x))
    xlabel(xlab, fontsize = 15)
    ylabel(ylab, fontsize = 15)
    xt = pyconvert(Vector{Int}, gca().get_xticks())
    xticks(xt, xt, fontsize = 12)
    yt = pyconvert(Vector{Int}, gca().get_yticks())
    yticks(yt, yt, fontsize = 12)
    xl = ax.get_xlim()
    yl = ax.get_ylim()
    r2 = round(1 - sum((y .- ỹ).^2) /sum((y .- mean(y)).^2), sigdigits = 3)
    α = round(lls[1], sigdigits = 2)
    β = round(lls[2], sigdigits = 2) 
    text(x = xl[0] + (xl[1] - xl[0])*0.03, y = yl[0] + (yl[1] - yl[0])*0.9, s = vary * " = " * string(α) * varx * " + " * string(β) * "\nR² = " * string(r2), fontweight = "bold", fontsize = 12)
    text(x = xl[0] + (xl[1] - xl[0])*0.85, y = yl[0] + (yl[1] - yl[0])*0.05, s = letter, fontweight = "bold", fontsize = 30)
end
tight_layout()
savefig(plotsdir("T_v_d18O.png"))


