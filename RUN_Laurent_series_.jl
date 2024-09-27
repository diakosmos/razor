include("Laurent.jl"); using .Laurent
using Plots
using LaTeXStrings

α=1; ν= 1/3; s=1
x_ = collect(range(-3,3,length=500))
σc_ = Laurent.σclassicalː.(x_;α=α,ν=ν)
fig = plot(x_,σc_,label="classical",linewidth=2)
xlabel!(L"x/L_v")
xlims!(-3,3)
display(fig)

# log10X = -2
f = ( 1/10^(3/8) )^(2/3)
α=1; ν= f* 1/3; s=f*1
#α=1; ν=1/3; s=100
xc = Laurent.xcː(α,ν,s)
xvm2 = Laurent.xvː(α,ν)
log10X = log10Xː(α,ν,s)
x_ = collect(range(-10*xc,10*xc,length=500))
σm2_ = Laurent.σː.(x_;α=α,ν=ν,s=s)
ϕm2_ = Laurent.ϕː.(x_;α=α,ν=ν,s=s)
fig = plot!(x_./xvm2,σm2_,label=L"\log_{10}X = -2",linestyle=:dash,linewidth=2)

# log10X = 0
f = ( 1/10^(3/8) )^(0/3)
α=1; ν= f* 1/3; s=f*1
#α=1; ν=1/3; s=100
xc = Laurent.xcː(α,ν,s)
xv0 = Laurent.xvː(α,ν)
log10X = log10Xː(α,ν,s)
x_ = collect(range(-10*xc,10*xc,length=500))
σ0_ = σː.(x_;α=α,ν=ν,s=s)
ϕ0_ = Laurent.ϕː.(x_;α=α,ν=ν,s=s)
fig = plot!(x_./xv0,σ0_,label=L"\log_{10}X = 0",linestyle=:dashdot,linewidth=2)


# log10X = 2
f = ( 1/10^(3/8) )^(-2/3)
α=1; ν= f* 1/3; s=f*1
#α=1; ν=1/3; s=100
xc = Laurent.xcː(α,ν,s)
xvp2 = Laurent.xvː(α,ν)
log10X = log10Xː(α,ν,s)
x_ = collect(range(-10*xc,10*xc,length=500))
σp2_ = σː.(x_;α=α,ν=ν,s=s)
ϕp2_ = Laurent.ϕː.(x_;α=α,ν=ν,s=s)
fig = plot!(x_./xvp2,σp2_,label=L"\log_{10}X = 2",linestyle=:dashdotdot,linewidth=2)
ylabel!(L"\Sigma / \Sigma_0")

display(fig)
savefig("Laurent.svg")
# To convert to eps (for example):
#"c:\Program Files\Inkscape\bin\inkscape.exe" Laurent.svg -o Laurent.eps --export-ignore-filters --export-ps-level=3
