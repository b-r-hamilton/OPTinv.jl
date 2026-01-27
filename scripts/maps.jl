
using PythonCall
if ! @isdefined(solutions)
    include("transientinversion.jl")
end

# ======== MAP OF SURFACE RECONSTRUCTION IN 70 year means ======= #
inds = collect(1100yr:200yr:2000yr)
stepsize = 70  
lev = -1.5:0.2:1.5

proj = ccrs.Orthographic(central_longitude=-80+55,central_latitude = 35)
noproj = ccrs.PlateCarree(central_longitude = 0)
fig, axs = subplots(nrows = length(inds), ncols = length(solutions),
                    subplot_kw = Dict("projection"=>proj), figsize = (8,12),
                    constrained_layout = true)
#(80,35,-80,30,5)
lon_box = [-80, 30]
lat_box = [35, 80]

for (i, sol) in enumerate(solutions) 
    Tu = sol isa solution ? Array(sol.ũ.dims[1]) : year.(lmrtime)yr
    Ty = sol isa solution ? Array(sol.y.dims[1]) : year.(lmrtime)yr 
    anom_index = findall(x->1850yr<x<1970yr, Tu)

    #remove 1850-1970anomaly from value we are plotting 
    if sol isa solution
        θ_ = [ustrip.(value.(x)) for x in sol.θ]
        θ = [x .- (sum(θ_[anom_index]) ./ length(anom_index)) for x in θ_]
    else
        θ = sol
        θ .-= (sum(θ[:, :, anom_index], dims = 3) ./ length(anom_index))
    end
            
    for (ii, y) in enumerate(inds)
        
        start_index = findall(x->x == y, Tu)[1]
        stop_index = findall(x->x == y + stepsize*yr, Tu)[1]
    
        if sol isa solution && y ∈ Ty
            plotme = γreshape(sum(θ[start_index:stop_index] ./ (stepsize/10)), sol.γ)' #100yr at 10yr res
        elseif sol isa Array
            plotme = sum(θ[:, :, start_index:stop_index], dims = 3)[:, :, 1] ./ stepsize #100yr at 1yr res
            plotme = plotme'
        end
        
        ax = axs[ii-1, i-1]
        if (sol isa solution && y ∈ Ty) || (sol isa Array)
            
            ax_hdl = ax.plot(orthographic_axes(lat_box..., lon_box..., 5)...,
                             color="black", linewidth=0.5,
                             transform=noproj)
            tx_path = first(ax_hdl)._get_transformed_path()
            
            path_in_data_coords, _ = tx_path.get_transformed_path_and_affine()
            polygon1s = mpath.Path(path_in_data_coords.vertices)
            ax.set_boundary(polygon1s)
            ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "110m", edgecolor="k", facecolor="gray"))

            ax.set_extent([lon_box..., lat_box[1]- 10, lat_box[2]], noproj)
            
            #add a cyclic point
            lon = sol isa solution ? ustrip.(sol.γ.lon) : lmrlon 
            #plotme, lon_ = cu.add_cyclic_point(plotme, coord = lon)
            lat = sol isa solution ? ustrip.(sol.γ.lat) : lmrlat
            lat_box_= [lat_box[1] - 20, lat_box[2] + 20]
            lon_box_ = [360-80 - 20, 30 + 20]
            plotme = geospatialsubset(plotme', lat, lon, lat_box_, lon_box_)
            plotme = plotme'
            cf = ax.contourf(lon, lat, plotme, cmap = "coolwarm", levels = lev, transform = noproj)
            c = ax.contour(lon, lat, plotme, colors = "black", levels = lev, transform = noproj)
            
            ax.clabel(c, fontsize = 8)
            gl = ax.gridlines(draw_labels = true)
            gl.top_labels = false
            gl.right_labels = false
            gl.bottom_labels = false
        else
            ax.clear()
            ax.axis("off") 
        end
    end
    tight_layout()
    savefig(plotsdir("surfacesol" * suffix * ".png"), bbox_inches = "tight",dpi = 600)
end



# ==================== SLOPE AT EVERY POINT, COMPARE TO OC2k SLOPES ==== #
sol = oldc
Tu = Array(sol.ũ.dims[1])
ind = findall(x->1900yr>x>1140yr, Array(Tu))

θmat = value.(ustrip.(hcat(sol.θ...)))
lls = Vector{Float64}(undef, size(θmat)[1])

for i in 1:size(lls)[1]
    lls[i] = linearleastsquares(ustrip.(Array(Tu)[ind]), θmat[i, ind])[1][1]
end
proj = ccrs.Orthographic(central_longitude=-80+55,central_latitude = 35)
noproj = ccrs.PlateCarree(central_longitude = 0)
fig, axs = subplots(nrows = 1, ncols = 1,
                    subplot_kw = Dict("projection"=>proj), figsize = (10,8))
ax = axs
ax_hdl = ax.plot(orthographic_axes(lat_box..., lon_box..., 5)...,
                 color="black", linewidth=0.5,
                 transform=noproj)
tx_path = ax_hdl[0]._get_transformed_path()
path_in_data_coords, _ = tx_path.get_transformed_path_and_affine()
polygon1s = mpath.Path(path_in_data_coords.vertices)
ax.set_boundary(polygon1s)
ax.add_feature(cfeature.NaturalEarthFeature("physical", "land", "110m", edgecolor="k", facecolor="gray"))
lon = ustrip.(oldc.γ.lon)
lat = ustrip.(oldc.γ.lat)
lat_box_= [lat_box[1] - 20, lat_box[2] + 20]
lon_box_ = [360-80 - 20, 30 + 20]
plotme = γreshape(lls .* 1000, oldc.γ)
@show NaNMath.maximum(plotme)
#plotme = geospatialsubset(plotme, lat, lon, lat_box_, lon_box_)
plotme, lon_ = cu.add_cyclic_point(plotme', coord = lon)

mpl = pyimport("matplotlib")

cf = ax.contourf(lon_, lat, plotme, cmap = "coolwarm", vmin = -3, vmax = 3, transform = noproj)
c = ax.contour(lon_, lat, plotme, colors = "black", vmin = -3, vmax = 3, transform = noproj) 
ax.clabel(c)

#find LLS of Ocean2K records in region of interest, if they have data between 1000 and 1500
lons = [-50, 20]
lats = [50, 90] 
oc2k_binned, oc2k_binnedv, oc2k_ages = loadOcean2kBinned()
oc2k, oc2kdata = loadOcean2k()

lon = oc2k[!, "lon"]
lat = oc2k[!, "lat"]
depth = oc2k[!, "depth"]
recinds = intersect(findall(x->-50 < x < 20, lon), findall(x->50<x<90, lat))
recs = oc2k[recinds, "name"]

function getrid(x)
    x[x.== "NaN"] .= NaN
    return convert(Vector{Float64}, x)
end

d = Dict() 
for (i, r) in enumerate(recs)
        
    locmatind = recinds[i]
    lat_ = lat[locmatind]
    lon_ = lon[locmatind]
    @show r
    inds = findall(x->occursin(r[begin:20], x), names(oc2kdata))
    cutoff = findall(x->x=='_', names(oc2kdata)[inds[1]])[1]
    lab = names(oc2kdata)[inds[1]][begin:cutoff-1]
    age = getrid(oc2kdata[!, inds[1]])
    sst = getrid(oc2kdata[!, inds[2]])
    post1000 = findall(x->1140<x<1900, age)
    #if NaNMath.maximum(age) > 1500
        lls, C = linearleastsquares(age[post1000], sst[post1000])
        @show lls[1] * 1000
        d[r] = (lat_, lon_, lls[1] * 1000)
    #end
end
oc2kms = 150
not_thornalley = [k for k in keys(d)][keys(d) .!= "Atlantic0220Thornalley2009__RAPiD-12-1K"]
xoc2k = Vector{Float64}(undef, length(not_thornalley))
yoc2k = Vector{Float64}(undef, length(not_thornalley))
coc2k= Vector{Float64}(undef, length(not_thornalley))


for (i, k) in enumerate(not_thornalley)

        xoc2k[i] = d[k][2]
        yoc2k[i] = d[k][1]
        coc2k[i] = d[k][3]
        

end
#[text(x = xoc2k[i] -4, y = yoc2k[i]-1, s = round.(coc2k, sigdigits = 1)[i], color = "white", zorder = 1000000000, transform = noproj, fontsize = 8) for i in 1:length(not_thornalley)]

s = scatter(xoc2k, yoc2k, c =coc2k , vmin = -3, vmax = 3, cmap = "coolwarm", s = oc2kms, edgecolors = "white", transform = noproj, zorder = 100000,linewidths = 2)
cb = colorbar(s, fraction = 0.025, pad = 0.0000001, orientation = "horizontal", location = "top")
cb.set_label("ΔT/Δt [K/kyr]", fontsize = 15)

key1 = "Arctic1147Bonnet2010_JM-06-WP-04-MCB"
key2 = "Arctic1148Spielhagen2011_MSM5/5-712"

MS = pyimport("matplotlib.markers").MarkerStyle

scatter(d[key1][2], d[key1][1], c = d[key1][3], vmin = -3, vmax = 3,  cmap = cm.balance, s = oc2kms, edgecolors = "white", transform = noproj, zorder = 100000,linewidths = 2, marker = MS("o", fillstyle = "right"))
scatter(d[key2][2], d[key2][1], c = d[key2][3], vmin = -3, vmax = 3,  cmap = cm.balance, s = oc2kms, edgecolors = "white", transform = noproj, zorder = 100000,linewidths = 2, marker = MS("o", fillstyle = "left"))

gl = ax.gridlines(draw_labels = true)
gl.top_labels = false
gl.right_labels = false
gl.bottom_labels = false

tight_layout()
savefig(plotsdir("ocean2kcomp.png"), dpi = 600)
#=
# ===== A SUPPLEMENTAL FIGURE ========== #
GS = pyimport("matplotlib.gridspec").GridSpec
subplotdir = plotsdir("oceank2compsub")
!isdir(subplotdir) && mkdir(subplotdir)
               
corenames = Dict(recs .=> ["MD95-2011", "ODP-162-984", "ENAM9606,M200309", "MD99-2275", "RAPiD-21-3K", "RAPiD-12-1K", "JM-06-WP-04-MCB", "MSM5-5-712"])

for (i, r) in enumerate(recs)
    figure()        
    inds = findall(x->occursin(r[begin:20], x), names(oc2kdata))
    cutoff = findall(x->x=='_', names(oc2kdata)[inds[1]])[1]
    lab = names(oc2kdata)[inds[1]][begin:cutoff-1]
    age = getrid(oc2kdata[!, inds[1]])
    sst = getrid(oc2kdata[!, inds[2]])
    post1000 = findall(x->1140<x<1900, age)
    lls, C = linearleastsquares(age[post1000], sst[post1000])
    plot(age, sst, label = r, color = "black", linewidth = 5)
    plot(age, age .* lls[1] .+ lls[2], color = "red", linewidth = 5) 
    xlim(1140,1900)
    xticks(1200:200:1800, 1200:200:1800, fontsize = 15)
    yt = gca().get_yticks()
    yticks(yt, yt, fontsize = 15)
    ytickrange = yt[length(yt)-1] - yt[0]
    text(x = 1200, y = yt[0] .+ 0.9ytickrange, s = string(round(lls[1] * 1000, sigdigits = 3)) * " °C/kyr", fontsize = 20)
    title(corenames[r], fontsize = 20)
    ylabel("Temperature [°C]", fontsize = 20)
    tight_layout()
    savefig(joinpath(subplotdir, corenames[r]), dpi = 600)
end



=#
