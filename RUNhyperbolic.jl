include("godunov.jl")
using .Godunov
using Plots
fig = Plots.plot()

"""
NN = 100 # no of cells
Nt = 64000000 # no of time steps
Np = 6400000 # how often to plot
Ns = 5000 # when to start plotting
"""
NN = 100 # no of cells
Nt = 8000 # no of time steps
Np = 1000 # how often to plot
Ns = 1#5000 # when to start plotting


r=Godunov.ring(NN,xi=0.5,xf=20.0)
Godunov.getAlphas!(r)
for i in 1:Nt
    if (mod(i,Np)==0) & (i>Ns)
        t = r.t
        Plots.plot!(r.X_[1:NN÷1],r.Σ_[1:NN÷1],label="t: $t",legend=:bottomright)#,title="Σ")#;legend=false)
        display(fig)
        Plots.plot!(r.X_[1:NN÷1],r.J_[1:NN÷1],label=false)#,label="t: $t",legend=:topleft)#,title="J")#;legend=false)
        display(fig)
    end
    Godunov.step!(r)
end

print("Done!")
