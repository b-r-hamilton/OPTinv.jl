#=
MODE MAGNITUDES THROUGH TIME 
=#
if ! @isdefined(solutions)
    include("transientinversion.jl")
end

# ====================== MODE 1 ====================================== #
gkw = ("height_ratios" => [1,1],)
fig, axs = subplots(2,1, gridspec_kw = gkw)
function specialplot(da::DimArray, ax, color::String; fb = true)
    x = ustrip.(Array(dims(da)[1]))
    y = Array(da)
    yunc = ustrip.(Measurements.uncertainty.(y))
    y = ustrip.(value.(y))
    ax.plot(x,y, color = color)
    if fb
        ax.fill_between(x, y1 = y .- yunc, y2 = y.+yunc, color = color, alpha = 0.5)
    end
    
end

specialplot(allc.u₀.x[:, At(1), At(:θ)], axs[0], "gray")
specialplot(allc.ũ.x[:, At(1), At(:θ)], axs[0], "red")
specialplot(oldc.ũ.x[:, At(1), At(:θ)], axs[0], "blue")
axs[0].set_ylabel("Mode Mag. [K]", fontsize = 15)
yt = axs[0].get_yticks()
axs[0].set_yticks(yt,yt,fontsize = 12)
xt = collect(500:500:2000)
axs[0].set_xticks(xt, fill("", length(xt)))
axs[0].set_xlabel("")
specialplot(allc.u₀.x[:, At(1), At(:θ)], axs[1], "gray", fb = false)
specialplot(allc.ũ.x[:, At(1), At(:θ)], axs[1], "red", fb = false)
specialplot(oldc.ũ.x[:, At(1), At(:θ)], axs[1], "blue", fb = false)
axs[1].set_ylabel("Mode Mag. [K]", fontsize = 15)
axs[0].text(x = 500, y = -60, s = "A", fontsize = 30, weight = "bold")
axs[1].text(x = 500, y = -20, s = "B", fontsize = 30, weight = "bold")
xticks(xt, xt, fontsize = 12)
axs[1].set_xlabel("Time [years CE]", fontsize = 15)
tight_layout()
savefig(plotsdir("u0utilde_mode1.png"))

# ==================== ALL MODES, WITH UNC. =========================== #
figure(figsize = (10,10))
for i in 1:11
    ax = subplot(4,3,i)
    plot(allc.u₀.x[:, At(i), At(:θ)], color = "gray")
    plot(allc.ũ.x[:, At(i), At(:θ)], color = "red")
    plot(oldc.ũ.x[:, At(i), At(:θ)], color = "blue")
    title(string(i), fontsize = 15)
    yt = ax.get_yticks()
    ax.set_yticks(yt,yt,fontsize = 12)
    if i ∈ [1,4,7,10]
        ylabel("Mode Mag. [K]", fontsize = 15)
    end
    xt = collect(500:500:2000)    
    if i ∉ [9,10,11]
        xticks(xt,fill("", length(xt)))
        xlabel("")
    else
        xlabel("Time [years CE]", fontsize = 15)

        xticks(xt, xt, fontsize = 12)
    end    
end
#tight_layout()
savefig(plotsdir("u0utilde.png"))


# ===================== FIRST THREE MODES, AS Z-SCORE, NO ERROR REP. ==== # 
for s in [:θ, :δ]
figure(figsize = (8,8))
    for (i, sol) in enumerate(solutions)
        subplot(length(solutions),1, i)
        modes = sol.ũ.dims[2]
        for m in 1:3
            #color = cmap.get_cmap("nipy_spectral")(1/(maximum(modes) - minimum(modes)) * m)
            color = ["purple", "green", "darkorange"][m]
            ls = ["solid", "dashed", "dotted"][m]
            tr = (minimum(sol.y.dims[1]), maximum(sol.y.dims[1]))
            ts = value.(sol.ũ.x[:, At(m), At(s)])
            plot((ts .- mean(ts)) ./ std(ts), label = string(m), color = color, zorder = length(modes) - m, linewidth = 4, linestyle = ls)
            
        end
        if i == 2

            xticks(600:200:2000, fontsize = 12)
            xlabel("Time [years CE]", fontsize = 15)
        else
            xticks(600:200:2000, fill("", length(600:200:2000)))
            xlabel("")
        end
        
        yticks(-4:2:6, fontsize = 12) 
        ylabel("Mode Magnitude Anomaly [σ]", fontsize = 15)
        #legend(loc = "center left")
        text(x = 485, y = 5, s = ["A", "B"][i], fontsize = 30, fontweight = "bold")
        grid()
        
        xlim(475, 2000)
        ylim(-4, 6)

        ax = gca()
        twin2 = ax.twinx()
        twin2.spines.right.set_position(("axes", 1.15))
        steinhilber = loadSteinhilber2009()
        st_time = 1950 .- steinhilber[!, "YearBP"]
        st = Measurements.measurement.(steinhilber[!, "dTSI"], steinhilber[!, "dTSI_sigma"])
        twin2.plot(st_time, steinhilber[!, "dTSI"], color = "grey")
        twin2.fill_between(x = st_time, y1 = steinhilber[!, "dTSI"] .-  steinhilber[!, "dTSI_sigma"], y2 = steinhilber[!, "dTSI"] .+  steinhilber[!, "dTSI_sigma"], color = "grey", alpha = 0.2)
        twin2.set_ylabel("Total Solar Insolation [Wm⁻²]", fontsize = 15, color = "grey")
        twin2.set_yticks(labels = -1.5:0.5:1.5, ticks = -1.5:0.5:1.5, fontsize = 12, color = "grey")
    end
    tight_layout()
    savefig(plotsdir("modemags"*string(s)* suffix* ".png"))
end
