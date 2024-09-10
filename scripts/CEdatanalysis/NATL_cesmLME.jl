#=
This script will calculate the mean SST for the following regions
- regional mean: 50W-20E, 50N-90N
- SPNA: 60W-0, 50-65N 
- Nordic Seas: 20W-20E, 65N-80N

For the following CESM-LME members
- each of the 13 ensemble members
- each of the 3 orbital runs 
- each of the 5 volcanically forced runs

Can use following command to move to casper
cd [path to OPTinv.jl/data]
scp NATL_cesmLME.jl [USERNAME]@casper.hpc.ucar.edu:/glade/u/home/[USERNAME]

Then log into casper, and run script with
`module load julia`
`julia NATL_cesmLME.jl`

Then, locally, run the following to pull the jld2 file down from casper
`scp [USERNAME]@casper.hpc.ucar.edu:/glade/u/home/[USERNAME]/NATL_cesmLME.jld2 .`

Can then run analysis script `NATL_cesmLME_analysis.jl
=#
using NCDatasets, Dates, JLD2, Statistics

ts = collect(DateTimeNoLeap(850, 2,1):Month(1):DateTimeNoLeap(2006, 1,1))

exp = vcat("00" .* string.(1:13),
           "ORBITAL.00" .* string.(1:3),
           "VOLC_GRA.00" .* string.(1:5))

sstdir = "/glade/campaign/cesm/collections/cesmLME/CESM-CAM5-LME/ocn/proc/tseries/monthly/SST"

#find the spatial indices of interest for each region 
templatefile = "b.e11.BLMTRC5CN.f19_g16.001.pop.h.SST.090001-099912.nc"
nc = NCDataset(joinpath(sstdir, templatefile))
lon = nc["TLONG"][:, :] #degrees E 
lat = nc["TLAT"][:, :] #degrees N

sst = nc["SST"][:, :, 1, :]; #2.5 secs
#we want to find which indices are missing 
notmissing = findall(x->!ismissing(x), sst[:, :, 1])
reg_ind = intersect(findall(x->360-50 < x || x < 20, lon),
                    findall(x->50<x<90, lat), notmissing)
spna_ind = intersect(findall(x->360-60 < x, lon),
                     findall(x->50<x<65, lat), notmissing)
ns_ind = intersect(findall(x->360-20 < x || x < 20, lon),
                   findall(x->65<x<80, lat), notmissing)
inds = [reg_ind, spna_ind, ns_ind]
mat = fill(NaN, (length(exp), length(inds), length(ts)))

close(nc)

files = readdir(sstdir) 

for (i, expid) in enumerate(exp)
    fileprefix = "b.e11.BLMTRC5CN.f19_g16."*expid*".pop.h.SST."
    #which files exist in this experiment id? 
    file_ind = findall(x->x ==1, occursin.(fileprefix, files))
    for fi in file_ind
	println("opening: " * files[fi])
        nc = NCDataset(joinpath(sstdir, files[fi]))
        sst = nc["SST"][:, :, 1, :]
        nctime = nc["time"][:]
        timeinds = findall(x->x∈ nctime, ts)
        for (j, spatialinds) in enumerate(inds)
            μ = [mean(sst[:, :, i][spatialinds]) for i in 1:size(sst)[3]]; #0.28s
            mat[i, j, timeinds] = μ
        end
        close(nc) 
    end
end

jldsave("/glade/u/home/hamiltonb/NATL_cesmLME.jld2"; mat, exp, inds, ts)  
