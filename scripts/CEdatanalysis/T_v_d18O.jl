#=
Analyzes TMI output for T v. S and T v. d18O from TMI output 
=#
import Pkg;Pkg.activate("../../")
using TMI, NCDatasets, PythonPlot, OPTinv, UnitfulLinearAlgebra, PythonCall, Statistics, DrWatson, PythonPlotExt, Measurements

TMIversion = "modern_180x90x33_GH11_GH12"
if ! @isdefined A
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
end
nc = NCDataset(TMIfile)
lat = γ.lat
lon = γ.lon 
xlabels = ["S [psu]", "δ¹⁸O [‰]", "δ¹⁸O [‰]"]
ylabels = ["T [°C]", "S [psu]", "T [°C]"]
vy = ["T", "S", "T"]
vx = ["S", "δ¹⁸O", "δ¹⁸O"]
letters = ["A", "B", "C"]
latpoi = 51.6
lonpoi = 360 - 21.9
lati = findmin(abs.(latpoi .- lat))[2]
loni = findmin(abs.(lonpoi .- lon))[2]

d = Dict()

d18O = nc["δ¹⁸Ow"][loni, lati, :] #.± nc["σδ¹⁸Ow"][loni, lati, :]
T = nc["θ"][loni, lati, :] #.± nc["σθ"][loni, lati, :]

S = nc["Sp"][loni, lati, :] #.± nc["σSp"][loni, lati, :]
pairs = [[S, T], [d18O, S]]


sdx = [3,2,2]
sdy = [2,3,2]
figure(figsize = (9, 4))

for (i, (p, xlab, ylab, varx, vary, letter)) in enumerate(zip(pairs, xlabels, ylabels, vx, vy, letters))
    ax = subplot(1,2,i)
    dind = (!).(isnan.(vec(p[1])))
    scatter(p[1][dind], p[2][dind], color = "black")

    
    x = vec(p[1])[dind]
    y = vec(p[2])[dind]
    lls, llsc = linearleastsquares(value.(x),value.(y))
    ỹ = value.(x).*lls[1] .+ lls[2]
    redind = findall(x->x>minimum(y), ỹ)
    plot(value.(x[redind]), ỹ[redind], color = "red", linewidth = 2)
    ỹ06 = (value.(x) .* 1/0.5)  .+ lls[2] .+ 0.17
    ỹ022 = (value.(x) .* 1/0.22)  .+ lls[2] .- 0.8
    @show length(x)
    if i == 2
        plot(value.(x[redind]), ỹ06[redind], color = "red", linewidth = 2, linestyle = "dashed")
        #plot(value.(x[redind])[10:end-5], ỹ022[redind][10:end-5], color = "red", linewidth = 2, linestyle = "dashdot")
    end
    
    r2 = round(1 - sum((y .- ỹ).^2) /sum((y .- mean(y)).^2), sigdigits = 3)
    α = round(lls[1], sigdigits = 2)
    β = round(lls[2], sigdigits = 2)


    xlabel(xlab, fontsize = 15)
    ylabel(ylab, fontsize = 15)
    #ylim(minimum(y), maximum(y))
    #xlim(minimum(x), maximum(x))
    
    xt = gca().get_xticks()
    xt = round.(pyconvert(Vector{Float64}, gca().get_xticks()), sigdigits = sdx[i])
    xticks(xt, xt, fontsize = 12)
    yt = round.(pyconvert(Vector{Float64}, gca().get_yticks()), sigdigits = sdy[i])
    yticks(yt, yt, fontsize = 12)
    xl = ax.get_xlim()
    yl = ax.get_ylim()
    text(x = xl[0] + (xl[1] - xl[0])*0.03, y = yl[0] + (yl[1] - yl[0])*0.9, s = vary * " = " * string(α) * varx * " + " * string(β) * "\nR² = " * string(r2), fontweight = "bold", fontsize = 12, color = "black")
    text(x = xl[0] + (xl[1] - xl[0])*0.85, y = yl[0] + (yl[1] - yl[0])*0.05, s = letter, fontweight = "bold", fontsize = 30)
    
end


tight_layout()
savefig(plotsdir("T_v_d18O.png"))
