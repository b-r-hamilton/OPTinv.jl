# OPTinv

[![Build Status](https://github.com/b-r-hamilton/OPTinv.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/b-r-hamilton/OPTinv.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/b-r-hamilton/OPTinv.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/b-r-hamilton/OPTinv.jl)
# What is this?
The source code for `OPTinv` holds methods to compute spatiotemporal surface maps of variables that explain the time evolution of some observed variable from benthic sediment core compilations. The scripts developed here focusing on applying this methodology to infer surface values of temperature and $\mathrm{\delta}^{18}\mathrm{O}\_{\mathrm{seawater}}$ from the [Lu et al, 2023](https://www.science.org/doi/10.1126/science.adf1646) dataset of 11 sediment core records of $\mathrm{\delta}^{18}\mathrm{O}\_{\mathrm{calcite}}$ from the Reykjanes Ridge in Iceland. 

# How to run 
Clone this repository, start a Julia session, activate the Julia environment in `OPTinv.jl`, and run `ex3.transientinversion.jl` to produce the two solutions detailed in the manuscript. Output can be plotted using any of the follow-on scripts. 

# Philosophy, and important dependencies 
The code developed here attempts to preserve information about the dimensionality and units of all variables. This is accomplished by using the following Julia libraries 
- [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl) to manage quantities, 
- [`DimensionalData`](https://github.com/rafaqz/DimensionalData.jl) for managing dimensioned data, 
- [`BLUEs.jl`](https://github.com/ggebbie/BLUEs.jl) for solving the linear underdetermined system presented here, 
- [`UnitfulLinearAlgebra.jl`](https://github.com/ggebbie/UnitfulLinearAlgebra.jl) for properly propagating physical quantities through linear algebra 

Information about the modern-day circulation is obtained from the Total Matrix Intercomparison (TMI) results from [Gebbie, 2010](https://journals.ametsoc.org/view/journals/phoc/40/8/2010jpo4272.1.xml) and [Gebbie, 2012](https://journals.ametsoc.org/view/journals/phoc/42/2/jpo-d-11-043.1.xml). These results are accessed through the following Julia libraries 
- [`TMI.jl`](https://github.com/ggebbie/TMI.jl)
- [`TMItransient.jl`](https://github.com/ggebbie/TMItransient.jl) 


# Scripts 
Transient Inversion of Lu et al, 2023 records 
- `ex1.offsets.jl`: stand-alone script that generates depth-profile of effective modern-day d18Oc and compares to the full range of CE variability in recorded sediment core d18Oc 
- `ex3.plotmodes.jl`: stand-alone script that plots surface spatial mode plots, impulse response of each mode at 3 select cores, and singular values w.r.t. sediment core/depth 
- `ex3.transientinversion.jl`: generates the 2 solutions detailed in the manuscript
- `ex3.transientinversionoc2k.jl`: attempt at constraining two solutions with Ocean2k cooling rate 
- `ex4.utilde.jl`: follow-on script, plots the mode solutions to `ex3` 
- `ex4.ytilde.jl`: follow-on script, plots the reconstructions at the core sites for each solution produced in `ex3`
- `ex5.regionmean.jl`: follow-on script, plots the regional means (for N. Atl. region and subregions) for solutions in `ex3`. Compares to Ocean2k and LMR. 
- `ex6.maps.jl`: follow-on script, plots the surface reconstruction maps for temperature 
- `ex7.effd18Oc.jl`: follow-on script, computes the effective d18Oc at the surface and compares it to the planktic stack from the EN539 sediment cores. 

Analysis of Common Era Data Products
- `CEdatanalysis/lmr.jl`: plots LMR and Mod-ERA maps
- `CEdatanalysis/NATL_cesmLME.jl`: a script to run on Casper that collects corresponding regional mean SSTs from relevant CESM-LME ensembles 
- `CEdatanalysis/NATL_cesmLME_analysis.jl`: script to analyze the mean SST timeseries collected in `NATL_cesmLME.jl` 
- `CEdatanalysis/ocean2k.jl`: plot all Ocean2k plots and compute the trend in Ocean2k values 
- `CEdatanalysis/en4MLD.jl`: Attempt at computing the bottom-of-MLD T and S in EN4 and comparing to my results. Probably needs better data quality and understanding of EN4 to make this any good 
- `CEdatanalysis/en4labrador.jl`: Hovmoller plot of Labrador Sea in EN4, T and S plots of Labrador Sea Water

Toy Problems 
- `toyproblems/ex2.steadystateinversion.jl`: `ex3`, but uses steady-state Gebbie 2010 assumption at each timestep. Hasn't been tested in a while, no idea if it works 
- `toyproblems/ex3.tests.jl`: toy problem 

Deglacial Problem 
- `deglac/lund2015.jl`

# Started package by doing
```
t = Template(; user = "b-r-hamilton", dir = "~/Code", authors = "B. Hamilton", julia = v"1.10", plugins = [License(; name = "MIT"), Git(; manifest = true, ssh=true), GitHubActions(; x86=false), Codecov(), Develop(), ],)
t("OPTinv.jl") 
```

