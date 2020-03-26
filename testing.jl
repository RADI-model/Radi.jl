module tst

using Colors, Plots
Plots.default(show=true)
Plots.closeall()

include("RADI.jl")
import .RADI

depths, oxy, poc = RADI.model(10.0, 1/128000, 128000)
cmap = colormap("RdBu", size(oxy)[2])

cs = 1
p1 = plot(depths*100, oxy[:, 1]*1e3, legend=false, c=cmap[cs])
p2 = plot(depths*100, poc[:, 1], legend=false, c=cmap[cs])
for sp in 2:size(oxy)[2]
    plot!(p1, depths*100, oxy[:, sp]*1e3, legend=false, c=cmap[cs+sp-1])
    plot!(p2, depths*100, poc[:, sp], legend=false, c=cmap[cs+sp-1])
end # for sp
plot(p1, p2, layout=(2, 1))

RADI.say_RADI()

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
