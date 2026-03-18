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
yt = -40:20:40
axs[0].set_yticks(yt,yt,fontsize = 12)
xt = collect(500:500:2000)
axs[0].set_xticks(xt, fill("", length(xt)))
axs[0].set_xlabel("")
specialplot(allc.u₀.x[:, At(1), At(:θ)], axs[1], "gray", fb = false)
specialplot(allc.ũ.x[:, At(1), At(:θ)], axs[1], "red", fb = false)
specialplot(oldc.ũ.x[:, At(1), At(:θ)], axs[1], "blue", fb = false)
axs[0].set_xlim(500,1980)
axs[1].set_ylabel("Mode Mag. [K]", fontsize = 15)
axs[0].text(x = 510, y = -40, s = "A", fontsize = 15, weight = "bold")
axs[1].text(x = 500, y = -17, s = "B", fontsize = 15, weight = "bold")
xlim(500,1980)
xticks(xt, xt, fontsize = 12)
axs[1].set_xlabel("Time [years CE]", fontsize = 15)
yt = -15:5:0
axs[1].set_yticks(yt,yt,fontsize = 12)
tight_layout()
savefig(plotsdir("u0utilde_mode1.png"), dpi = 600)

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
    else
        ylabel("")
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
tight_layout()
savefig(plotsdir("u0utilde.png"), dpi = 600)

