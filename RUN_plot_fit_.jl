include("Laurent.jl");  using .Laurent
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

x8_ = x8a_ / 0.033
fig = plot(x8_,Σ8a_,label="N-body",linewidth=2)

α = 1.0 #0.0316
ν = 1.0 /3
s = 3.1e1
x_ = range(-10,10,length=498)
xc = xcː(α,ν,s)
log10X = log10Xː(α,ν,s)
σ_ = σː.(x_;α=α,ν=ν,s=s)
xlims!(-4,4)
xlabel!(L"x/L_c")
ylabel!(L"Σ/Σ_0")
fig = plot!(x_./xc,σ_,linestyle=:dash,label=L"\log_{10}X = 4.5",linewidth=2)  #"*"$log10X")
#=
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
=#
display(fig)
savefig("Plotfit.pdf")