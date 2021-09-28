include("godunov_new.jl")
using .Godunov
using Plots
fig = Plots.plot()

# Physical parameters:
α = 1.0
ν = 1.0
s = 0.001

# Other params:
Nc = 40 # no of cells
tmax = 1.0e2 # max time
dt = 3.0e1 # how often to plot
tp = 0.0e0#5000 # when to start plotting

r=Godunov.ring(Nc, (α,ν,s))
Godunov.step!(r)

xh = r.Ls["Lcrit"]
x_ = (r.X_[1:end-1]+r.X_[2:end])/2.0
x2_ = 2*r.X_[end].-reverse(x_)
xx_ = cat(x_,x2_,dims=1)

#for i in 1:10
for tt in 0.0:dt:tmax
    tx = r.t; X_ = r.X_; Σ_ = r.Σ_; Φ_ = r.Φ_
    if tx ≥ tp
        ΣΣ_ = cat(Σ_,reverse(Σ_),dims=1)
        ΦΦ_ = cat(Φ_,-reverse(Φ_),dims=1)
        Plots.plot!(xx_,ΣΣ_,label="t: $tx",legend=false)#,xlim=(0,3))#:bottomright)#,title="Σ")#;legend=false)
        display(fig)
        print("time: $tx\n")
#        Plots.plot!(xx_,ΦΦ_,label=false,linestyle=:dash)#,label="t: $t",legend=:topleft)#,title="J")#;legend=false)
#        display(fig)
    end
#    for j in 1:800
        Godunov.step!(r)#,dt)
#    end
    Godunov.tstep!(r,dt)
end

tf=r.t
print("Done at t= $tf\n")
