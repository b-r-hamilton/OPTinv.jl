#=
RECONSTRUCTION OF SOLUTIONS AT CORE SITES 
=# 
if ! @isdefined(solutions)
    include("ex3.transientinversion_abstract.jl")
end

locs = core_locations() 
min_age = [minimum(s.y.dims[1]) for s in solutions]
max_age = [maximum(s.y.dims[1]) for s in solutions]

xl = [minimum(min_age), maximum(max_age)]

#what does the "old" solution look like, propagated to all core sites? 
ỹold = allc.predict(value.(oldc.ũ.x))
Cññold =allc.E*oldc.ũ.C*allc.E'
ỹold = DimEstimate(allc.E*oldc.ũ.v, parent(Cññold) .* permil^2, allc.y.dims)
n = UnitfulMatrix(allc.y.v) - ỹold.v
Cnn = UnitfulMatrix(ustrip.(allc.y.C), fill(permil, length(n)), fill(permil^-1, length(n)))
Jdata = n' * inv(Cnn) * n
rmserror = n' * n

fig = figure(figsize = (8, 10))
for (i, core) in enumerate(allcores) 
    subplot(4,3, i)
    #we want to plot y from the solution that has the longest record 
    @show core_in_sol = [core ∈ s.y.dims[2] for s in solutions]
    ma_mod = copy(min_age)
    ma_mod[(!).(core_in_sol)] .= 1980yr
    @show ma_mod
    @show plotfrom = findmin(ma_mod)[2]
    for (j, s) in enumerate(solutions)
        if j == plotfrom
            plot(s.y.x[:, At(core)], label = "y", color = "black", zorder = j)
            #plot(s.ỹ₀.x[:, At(core)], label = "ỹ₀", color = "green", zorder = j)
        end
        if core in s.y.dims[2]
            plot(s.ỹ.x[:, At(core)], label = "ỹ: " * s.name, color = s.color, zorder = j)
        end
        if core ∉ s.y.dims[2]
            plot(ỹold.x[:, At(core)], color = "blue", zorder = j, linestyle = "dashed")
        end 
    end
    ylim(-0.2, 0.2)
    xlim(ustrip.(xl))


    if i ∈ [9,10,11]
        xticks(1000:250:1750, fontsize = 12)
        xlabel("Time [years CE]", fontsize = 15)
    else
        xticks(1000:250:1750, ["", "", "", ""])
        xlabel("")
    end
    
    if i ∈ [1,4,7,10]
        yticks(-0.2:0.2:0.2, fontsize = 12) 
        ylabel(L"\mathrm{\delta}^{18}\mathrm{O}_\mathrm{calcite}" * " [‰]", fontsize = 15)
    else
        ylabel("")
        yticks(-0.2:0.2:0.2, ["", "", ""])
    end
        
    title(string(core) * ", " * string(locs[core][3]) * "m", fontsize = 15)
    #grid(zorder = 1)
    gca().invert_yaxis()
end
tight_layout()
savefig(plotsdir("reconsol" * suffix * ".png"))

