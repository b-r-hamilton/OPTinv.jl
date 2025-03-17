import Pkg;Pkg.activate("../../")
using TMI, NCDatasets, PythonPlot, OPTinv, UnitfulLinearAlgebra, PythonCall, Statistics, DrWatson, PythonPlotExt, Measurements

TMIversion = "modern_180x90x33_GH11_GH12"
if ! @isdefined A
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
end
nc = NCDataset(TMIfile)
lat = γ.lat
lon = γ.lon 
#=
natl = TMI.surfaceregion(TMIversion, "NATL").tracer
lab = TMI.surfaceregion(TMIversion, "LAB").tracer
gin = TMI.surfaceregion(TMIversion, "GIN").tracer


latgrid = repeat(lat', inner = (180, 1))
longrid = repeat(lon, outer = (1, 90))
Sp = nc["Sp"][:, :, 1]
θ = nc["θ"][:, :, 1]
d18O = nc["δ¹⁸Ow"][:, :, 1]
indsp = findall(x->x >34, Sp)
ind = (LAB = intersect(indsp, findall(x->x == 1, lab)), GIN = intersect(indsp, findall(x->x == 1, gin)))

figure(figsize = (10,10));
lon_box = [-80, 30]
lat_box = [60, 90]
proj = ccrs.PlateCarree()
noproj = ccrs.PlateCarree()
l = 15
vmins = [-1.5, 31.5, -2]
#vmaxs = [21, 38, 1.5]
step = [1, 0.5, 0.25]
sd = [2,3,3]

ind_core_lab = intersect(union(findall(x->2 < x < 5, θ),findall(x->34.9 < x < 35.05, Sp)), findall(x->x > -3, d18O), findall(x->x == 1, lab), findall(x->x >45, latgrid), findall(x-> 290 < x > 270, longrid))
ind_core_gin = intersect(union(findall(x->2 < x < 5, θ),findall(x->34.9 < x < 35.05, Sp)), findall(x->x > -3, d18O), findall(x->x == 1, gin),findall(x->x >45, latgrid))

for (i, mat) in enumerate([θ, Sp, d18O])
    ax = subplot(3,1,i;projection = proj);
    #ax.coastlines()
    lev = range(start = vmins[i], step = step[i], length = l)
    cf = contourf(γ.lon, γ.lat, mat', levels = lev, transform = noproj, alpha = 0.5);
    colorbar(cf)
    c = contour(γ.lon, γ.lat, mat', levels = lev, colors = "black", transform = noproj)
    scatter(γ.lon[[i[1] for i in ind_core_lab]], γ.lat[[i[2] for i in ind_core_lab]], color = "red", s = 1)
    scatter(γ.lon[[i[1] for i in ind_core_gin]], γ.lat[[i[2] for i in ind_core_gin]], color = "blue", s = 1)
    ax.clabel(c)
    ax.set_extent([lon_box..., lat_box[1]-10, lat_box[2]], noproj)
end



figure()
subplot(;projection = ccrs.PlateCarree())
scatter(γ.lon[[i[1] for i in ind_core_lab]], γ.lat[[i[2] for i in ind_core_lab]], color = "red")
scatter(γ.lon[[i[1] for i in ind_core_gin]], γ.lat[[i[2] for i in ind_core_gin]], color = "blue")
gca().coastlines()

figure()
scatter(d18O[ind_core_lab], θ[ind_core_lab], color = "red", label = "LAB")
title("Range of θ, Sp for sed. cores: 2 < θ < 5, 34.9 < Sp < 35.05")
scatter(d18O[ind_core_gin], θ[ind_core_gin], color = "blue", label = "GIN")
legend()
xlabel("d18o [‰]")
ylabel("θ [°C]")
=#
#=
indx = [γ.lon[i[1]] for i in ind]
indy = [γ.lat[i[2]] for i in ind]
scatter(indx, indy, color = "gray", marker = "x", s = 5)
=#
#=
ind = (LAB = ind_core_lab, GIN = ind_core_gin)

d = Dict()
for k in keys(ind)
    d18O = nc["δ¹⁸Ow"][:, :, 1][ind[k]]
    σd18O = nc["σδ¹⁸Ow"][:, :, 1][ind[k]]
    T = nc["θ"][:, :, 1][ind[k]]
    S = nc["Sp"][:, :, 1][ind[k]]

    pairs = [[S, T], [d18O, S], [d18O, T]]
    d[k] = pairs
end
=#
xlabels = ["S [psu]", "δ¹⁸O [‰]", "δ¹⁸O [‰]"]
ylabels = ["T [°C]", "S [psu]", "T [°C]"]
vy = ["T", "S", "T"]
vx = ["S", "δ¹⁸O", "δ¹⁸O"]
letters = ["A", "B", "C"]
#=
colors = (LAB = "red", GIN = "blue")
figure(figsize = (13, 4))
for k in keys(d)
    pairs = d[k]
    color = colors[k] 
    for (i, (p, xlab, ylab, varx, vary, letter)) in enumerate(zip(pairs, xlabels, ylabels, vx, vy, letters))
        ax = subplot(1,3,i)
        scatter(p[1], p[2], color = color, alpha = 0.4)
        #dind = (!).(isnan.(vec(p[1])))
        #=
        x = vec(p[1])#[dind]
        y = vec(p[2])#[dind]
        
        lls, llsc = linearleastsquares(x,y)
        ỹ = x.*lls[1] .+ lls[2]
        redind = findall(x->x>minimum(y), ỹ) 
        plot(x[redind], ỹ[redind], color = color, linewidth = 2)
        r2 = round(1 - sum((y .- ỹ).^2) /sum((y .- mean(y)).^2), sigdigits = 3)
        α = round(lls[1], sigdigits = 2)
        β = round(lls[2], sigdigits = 2)
        xl = ax.get_xlim()
        yl = ax.get_ylim()
        text(x = xl[0] + (xl[1] - xl[0])*0.03, y = yl[0] + (yl[1] - yl[0])*0.9, s = vary * " = " * string(α) * varx * " + " * string(β) * "\nR² = " * string(r2), fontweight = "bold", fontsize = 12, color = color)
        #text(x = xl[0] + (xl[1] - xl[0])*0.85, y = yl[0] + (yl[1] - yl[0])*0.05, s = letter, fontweight = "bold", fontsize = 30)
        =#
        xlabel(xlab, fontsize = 15)
        ylabel(ylab, fontsize = 15)
        #uylim(minimum(y), maximum(y))
#xlim(minimum(x), maximum(x))
    end
end
        xt = round.(pyconvert(Vector{Float64}, gca().get_xticks()), sigdigits = 3)
        xticks(xt, xt, fontsize = 12)
        yt = round.(pyconvert(Vector{Float64}, gca().get_yticks()), sigdigits = 3)
        yticks(yt, yt, fontsize = 12)

tight_layout()
savefig(plotsdir("T_v_d18O.png"))
=#
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
