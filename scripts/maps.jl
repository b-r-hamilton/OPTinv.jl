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
            cf = ax.contourf(lon, lat, plotme, cmap = cm.balance, levels = lev, transform = noproj)
            c = ax.contour(lon, lat, plotme, colors = "black", levels = lev, transform = noproj)
            ax.clabel(c)
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
cf = ax.contourf(lon_, lat, plotme, cmap = cm.balance, vmin = -3, vmax = 3, transform = noproj)
c = ax.contour(lon_, lat, plotme, colors = "black", vmin = -3, vmax = 3, transform = noproj) 
ax.clabel(c)

#find LLS of Ocean2K records in region of interest, if they have data between 1000 and 1500
lons = [-50, 20]
lats = [50, 90] 
oc2k_binned, oc2k_binnedv, oc2k_ages = loadOcean2kBinned()
oc2k, oc2kdata = loadOcean2k()

lon = oc2k[!, "lon"]
lat = oc2k[!, "lat"]
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
    if NaNMath.maximum(age) > 1500
        lls, C = linearleastsquares(age[post1000], sst[post1000])
        @show lls[1] * 1000
        d[r] = (lat_, lon_, lls[1] * 1000)
    end
end
for (i, k) in enumerate(keys(d))
    if d[k][1] < 70 #don't include disagreeing Arctic recs
        s = scatter(d[k][2], d[k][1], c = d[k][3], vmin = -3, vmax = 3, cmap = cm.balance, s = 100, edgecolors = "white", transform = noproj, zorder = 100000,linewidths = 2)
    end
end
gl = ax.gridlines(draw_labels = true)
gl.top_labels = false
gl.right_labels = false
gl.bottom_labels = false

tight_layout()
savefig(plotsdir("ocean2kcomp.png"), dpi = 600)

