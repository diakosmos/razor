include("godunov_3.jl");   using .Godunov 
using Plots 
using LaTeXStrings 
using CSV
using DataFrames

df8a = CSV.read("fig_8a_.csv", DataFrame; header=["x","Σ"], types=Float64)
x8a_ = df8a[!,"x"]
Σ8a_ = df8a[!,"Σ"]

fig = plot(x8a_,Σ8a_)
display(fig)

#################################

α = 1.0 #0.0316
ν = 1.0 /3
s = 1.0e1
styles_ = [:dot,:dashdot,:dash,:dashdotdot]

Nc = 180 # no of cells
tmax = 70.0e0 # max time
r=Godunov.ring(Nc, (α,ν,s), xfr=6.0, xexp=1)
Godunov.step!(r)

xh = r.Ls["Lcrit"]
xv = r.Ls["Lvisc"]
log10X = round(log10(r.Ls["X"]),sigdigits=3)
x_ = (r.X●[1:end-1]+r.X●[2:end])/2.0
x2_ = 0.0.-reverse(x_)
xx_ = cat(x2_,x_,dims=1)
Godunov.tstep!(r,tmax)

ΣΣ_ = cat(reverse(r.Σ_[2:end-1]),r.Σ_[2:end-1],dims=1)
fig = plot()
factor = 34
plot!(xx_./xh/factor,ΣΣ_/r.Σ_[end],xlim=(-0.15,0.15),ylim=(0,3.0), label=L"\log_{10} X = "*"$log10X",legend=:topright,linestyle=styles_[1],linewidth=1.0)#,title="Σ")#;legend=false)        display(fig)
plot!(x8a_,Σ8a_;style=styles_[2],label="N-body")
ylabel!(L"\Sigma / \Sigma_0")
xlabel!(L"x")
display(fig)