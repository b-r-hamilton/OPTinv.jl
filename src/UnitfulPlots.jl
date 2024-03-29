import PyPlot.plot
import PyPlot.scatter

function plot(y::DimArray; kwargs...)
    unc = ustrip.(uncertainty.(vec(y)))    
    x = ustrip.(collect(first(dims(y))))
    xunit = unit(first(first(dims(y))))
    yunit = unit(first(y))
    y = ustrip.(value.(vec(y)))
    f = fill_between(x = x, y1 = y .- unc, y2 = y .+ unc, alpha = 0.5; kwargs...)
    plot(x, y; kwargs...)
    xlabel(string(xunit))
    ylabel(string(yunit))
    return x, y 
end

function scatter(x::Vector{Quantity{Measurement{T}}}, y::Vector{Quantity{Measurement{T}}}; kwargs...) where T <: AbstractFloat
    xerr = uncertainty.(ustrip.(x))
    yerr = uncertainty.(ustrip.(y))
    xunit = unit(first(x))
    yunit = unit(first(y))
    x = value.(ustrip.(x))
    y = value.(ustrip.(y))
    errorbar(x, y, xerr = xerr, yerr = yerr; kwargs...)
    xlabel(string(xunit))
    ylabel(string(yunit))
end

