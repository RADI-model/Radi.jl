module tst

using Colors, Plots, Profile, ProfileView
Plots.default(show=true)
Plots.closeall()

include("RADI.jl")
import .RADI

stoptime = 5/8760 # 5.0
interval = 1/8760 # 0.5/128000
saveperXsteps = 1 # 2*128000
oxy_i = 0.0
poc_i = 0.0

function radiplot(oxy_i, poc_i)
    @time depths, oxy, poc = RADI.model(stoptime, interval, saveperXsteps,
        oxy_i, poc_i)
    ntps = size(oxy)[2]
    cmap = colormap("RdBu", ntps)
    cs = ntps
    p1 = plot(depths*100, oxy[:, 1]*1e3, legend=false, c=cmap[cs])
    p2 = plot(depths*100, poc[:, 1], legend=false, c=cmap[cs])
    for sp in 2:ntps
        plot!(p1, depths*100, oxy[:, sp]*1e3, legend=false, c=cmap[cs-sp+1])
        plot!(p2, depths*100, poc[:, sp], legend=false, c=cmap[cs-sp+1])
    end # for sp
    plot(p1, p2, layout=(2, 1))
    return depths, oxy, poc
end # function radiplot

showprofile = false
if showprofile
    Profile.clear()
    @profile depths, oxy, poc = RADI.model(stoptime, interval, saveperXsteps,
        oxy_i, poc_i)
    ProfileView.view()
    RADI.say_RADI()
else
    depths, oxy, poc = radiplot(oxy_i, poc_i);
end

# using Base.SimdLoop
#
# function forsimd(xmax::Int)
#     y = 0
#     @simd for x in 1:xmax
#         y += sqrt(x)
#     end # for x
#     return y
# end
#
# function fornormal(xmax::Int)
#     y = 0
#     for x in 1:xmax
#         y += sqrt(x)
#     end # for x
#     return y
# end
#
# xmax = 10000000
# @time forsimd(xmax)
# @time fornormal(xmax)

# abstract type Solid end
#
# getx(value) = value::Solid


# # Trying out composite type but slower and too clever for its own good
# struct PorewaterVariable
#     depths::Array{Float64,1}
#     start::Array{Float64,1}
#     previous::Array{Float64,1}
#     now::Array{Float64,1}
#     npts::Int64
#     overlying::Float64
# end # struct PorewaterVariable
#
# function PorewaterVariable(depths::Array{Float64,1}, start::Array{Float64,1},
#         overlying::Float64)
#     npts = length(depths)
#     previous = fill(NaN, npts)
#     PorewaterVariable(depths, start, previous, start, npts, overlying)
# end # function PorewaterVariable
#
# function PorewaterVariable(depths::Array{Float64,1}, start::Float64,
#         overlying::Float64)
#     npts = length(depths)
#     start = fill(start, npts)
#     previous = fill(NaN, npts)
#     PorewaterVariable(depths, start, previous, start, npts, overlying)
# end # function PorewaterVariable
#
# function now2previous!(myh::PorewaterVariable)
#     myh.previous[:] = myh.now
# end # function now2previous!

end # module tst
