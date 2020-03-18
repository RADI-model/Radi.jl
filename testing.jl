module Tst

using Plots
Plots.default(show=true)
Plots.closeall()

include("RADI.jl")
import .RADI

depths, oxy = RADI.model(0.0, 10.0, 0.1, 10)

plot(depths, oxy)

RADI.say_RADI()

end # module Tst
