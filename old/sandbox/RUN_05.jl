include("sand_05.jl")
using .Sand_05
using Plots
fig = Plots.plot()


NN = 160 # no of cells
tmax = 6.0 # max time
dt = 2.0 # how often to plot
tp = 0.0#5000 # when to start plotting


r= Sand_05.ring(NN; U=-2.0)#5.5)#5.5)
for t in 0.0:dt:tmax
#for i in 1:3#10
    t = r.t
    if t ≥ tp
        x0_ = 0.5 * (r.X_[1:end-1] + r.X_[2:end])
        p1 = Plots.plot!(x0_[1:end],r.Σ_[1:end],label="t: $t",legend=false)#:bottomright)#,title="Σ")#;legend=false)
        #Plots.scatter!(x0_[1:end],r.Σ_[1:end],label=false)
        display(fig)
        p2 = Plots.plot!(x0_[1:end],20*r.J_[1:end],label=false,linestyle=:dash)#,label="t: $t",legend=:topleft)#,title="J")#;legend=false)
        #Plots.scatter!(x0_[1:end],10*r.J_[1:end],shape=:+,label=false)
        display(fig)
    end
    Sand_05.tstep!(r,dt)
    #Sand_05.step!(r)
end
