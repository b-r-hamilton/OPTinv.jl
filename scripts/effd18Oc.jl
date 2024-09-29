#=
Compute the effective d18Oc at the surface 
=# 
if ! @isdefined(solutions)
    include("ex3.transientinversion_abstract.jl")
end

## ============ EFFECTIVE δ¹⁸Oc ============= ##
modes = allc.spatialmodes
cdims = vec(covariancedims(allc.ũ.dims))
Tu = Array(allc.ũ.dims[1])
Vsub = UnitfulMatrix(modes[:, regionmeanindices], fill(NoUnits, 11), fill(NoUnits, 500))
Ey = fill(1/500, 500)
Edag = Vsub * Ey

mat = Matrix{Float64}(undef,length(Tu),length(cdims))
mat .= 0
for (i,t) in enumerate(Tu)
    xindsθ = findall(x->x[1]==t && x[3] == :θ,cdims)
    xindsδ = findall(x->x[1]==t && x[3] == :δ,cdims)
    mat[i, xindsθ] = Edag * -0.224
    mat[i, xindsδ] = Edag
end

mat = UnitfulMatrix(mat, fill(permil, length(Tu)), unitrange(allc.ũ.v))
figure()
for (tit, sol) in zip(["old", "all"], solutions)
    Css = mat*sol.ũ.C*transpose(mat)
    y = UnitfulMatrix(vec(mat*sol.ũ.v))
    unc = UnitfulMatrix(sqrt.(diag(parent(Css))), unitrange(y))
    plot(DimArray(measurement.(vec(y), vec(unc)), Ti(Array(sol.ũ.dims[1]))), label = tit, color = sol.color)
    if tit == "old"
        ind = findall(x->1900yr>x>1140yr, Tu)
        display(parent(y)[ind])
        l, c= linearleastsquares(ustrip.(Array(Tu))[ind], parent(y)[ind], C = parent(Css)[ind,ind])
        @show l
        display(sqrt.(diag(c)))
    end
end

xlabel("Time [years CE]", fontsize = 15)
ylabel("Effective " * L"\mathrm{\delta}^{18}\mathrm{O_{calcite}}" * " [‰]", fontsize = 15)
xticks(600:200:2000, fontsize = 12)
yticks(-0.2:0.1:0.2, fontsize = 12)
gca().invert_yaxis()
tight_layout()
savefig(plotsdir("effd18Omean.png"))

# ============ EFFECTIVE d18Oc compared to PLANKTIC STACK ================= #
cdims = vec(covariancedims(allc.ũ.dims))
Tu = Array(allc.ũ.dims[1])
locs = core_locations()
corelons = 360 .+ [c[1] for c in locs]
corelats = [c[2] for c in locs]
surfind = γbox(oldc.γ, minimum(corelats)-1, maximum(corelats)+1, minimum(corelons)-1, maximum(corelons)+1)

Vsub = UnitfulMatrix(modes[:, surfind], fill(NoUnits, 11), fill(NoUnits, 2))
Ey = fill(1/2, 2)
Edag = Vsub * Ey

mat = Matrix{Float64}(undef,length(Tu),length(cdims))
mat .= 0
for (i,t) in enumerate(Tu)
    xindsθ = findall(x->x[1]==t && x[3] == :θ,cdims)
    xindsδ = findall(x->x[1]==t && x[3] == :δ,cdims)
    mat[i, xindsθ] = Edag * -0.224
    mat[i, xindsδ] = Edag
end

mat = UnitfulMatrix(mat, fill(permil, length(Tu)), unitrange(allc.ũ.v))
figure()
for (tit, sol) in zip(["old", "all"], solutions)
    Css = mat*sol.ũ.C*transpose(mat)
    y = UnitfulMatrix(vec(mat*sol.ũ.v))
    unc = UnitfulMatrix(sqrt.(diag(parent(Css))), unitrange(y))
    ty = sol.y.dims[1]
    plot(DimArray(measurement.(vec(y), vec(unc)), Ti(Array(sol.ũ.dims[1])))[ty[begin]..ty[end-1],:], label = tit, color = sol.color)
    if tit == "old"
        ind = findall(x->1800yr>x>1200yr, Tu)
        l, c= linearleastsquares(ustrip.(Tu)[ind], parent(y)[ind], C = parent(Css[ind,ind]))
        @show l
        display(sqrt.(diag(c)))
    end
end

#make a planktic stack
Gb = loadcores(oldcores, dir = OPTinv.datadir("runBacon"), rules = ["G.b"])
Gi = loadcores(oldcores, dir = OPTinv.datadir("runBacon"), rules = ["G.i"])
Gbunc = reshape(sqrt.(diag(parent(Gb.Cnn.mat)))permil, size(Gb.y))
Giunc = reshape(sqrt.(diag(parent(Gi.Cnn.mat)))permil, size(Gi.y))
mat = hcat(Matrix(Gi.y) .± Giunc, Matrix(Gb.y) .± Gbunc)
mat = Matrix(Gb.y) .± Gbunc
inds = findall(x->1980yr > x > 1850yr, Array(Gb.y.dims[1]))
#remove 1850-1980 mean
mat .-= mean(value.(mat[inds, :]), dims = 1)
plot(DimArray(mean(mat, dims = 2)[:], Ti(Array(Gb.y.dims[1]))), color = "black")
xlabel("Time [years CE]", fontsize = 15)
ylabel("Effective " * L"\mathrm{\delta}^{18}\mathrm{O_{calcite}}" * " [‰]", fontsize = 15)
xticks(600:200:2000, fontsize = 12)
yticks(-0.4:0.1:0.3, fontsize = 12)
gca().invert_yaxis()
tight_layout()
savefig(plotsdir("effd18Ostack.png"))
