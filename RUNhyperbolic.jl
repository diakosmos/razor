include("godunov.jl")
using .Godunov
using Plots
fig = Plots.plot()
r=Godunov.ring()
Godunov.getAlphas!(r)
for i in 1:100
    Plots.plot!(r.X_,r.Σ_;legend=false)
    display(fig)
    Godunov.step!(r)
end
