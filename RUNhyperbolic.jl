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
NN = 40 # no of cells
#Nt = 8000 # no of time steps
#Np = 1000 # how often to plot
#Ns = 1#5000 # when to start plotting
tmax = 1.0e3 # max time
dt = 1.0e2 # how often to plot
tp = 0.0e0#5000 # when to start plotting


r=Godunov.ring(NN)#,xi=0.4,xf=3.0)
Godunov.getAlphas!(r)
#for i in 1:10
for tt in 0.0:dt:tmax
    tx = r.t
    if tx ≥ tp
        Plots.plot!(r.X_[1:NN÷1],r.Σ_[1:NN÷1],label="t: $tx",legend=false,xlim=(0,3))#:bottomright)#,title="Σ")#;legend=false)
        display(fig)
        print("time: $tx\n")
        Plots.plot!(r.X_[1:NN÷1],r.Φ_[1:NN÷1],label=false,linestyle=:dash)#,label="t: $t",legend=:topleft)#,title="J")#;legend=false)
        display(fig)
    end
#    for j in 1:800
        Godunov.step!(r)#,dt)
#    end
    Godunov.tstep!(r,dt)
end
tf=r.t
print("Done at t= $tf\n")
