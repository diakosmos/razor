include("classical.jl")
using .Classical
using Plots

r = Classical.ring(32;xf=2.0)
fig = Plots.plot()

for b in [0.0,1.0,2.0]#,-0.5,1.0,2.0,3.0]
    r.β = b
    for j in 1:1000000
        Classical.iterate!(r)
    end
    Plots.plot!(r.x_[1:32],r.σ_[1:32];legend=false)
    display(fig)
end
