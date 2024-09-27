include("godunov.jl");   using .Godunov
#include("godunov_full_.jl");   using .GodunovFull 
Godunov = Godunov
using Interpolations
using Plots
using LaTeXStrings

PLOT = true

# Physical parameters:
α = 1.0 #316.0 #0.0316
ν = 1.0e0 /3 #10.0 /3
s_ = (1.0e-2,) #(1.0e-2,1.0e-1,1.0,1.0e1)
styles_ = [:dot,:dashdot,:dash,:dashdotdot]

#title!("α: $α   ν: $ν    s: $s")

# Other params:
Nc = 20 # no of cells
tmax = 1.0e2 # max time
dt = 3.0e1 # how often to plot
tp = 0.0 #e2#5000 # when to start plotting

for s in s_
    if PLOT
        fig = Plots.plot()
        xlabel!(L"x / L_v")
        ylabel!(L"\Sigma/\Sigma_0")
    end
    #style = pop!(styles_)
    r=Godunov.ring(Nc, (α,ν,s), xf=12.0, xexp=1)
    Godunov.step!(r)

    xh = r.Ls["Lcrit"]
    xv = r.Ls["Lvisc"]
    log10X = round(log10(r.parmd["X"]),sigdigits=3)
    x_ = (r.X●[1:end-1]+r.X●[2:end])/2.0
    #x2_ = 2*r.X●[end].-reverse(x_)
    x2_ = 0.0.-reverse(x_)
    #xx_ = cat(x_,x2_,dims=1)
    xx_ = cat(x2_,x_,dims=1)
    Godunov.step!(r)
    #for i in 1:10
    for tt in 0.0:dt:tmax
        tx = r.t; Σ_ = r.Σ_; Φ_ = r.Φ_
        #if tx ≥ tp
            #ΣΣ_ = cat(Σ_[2:end-1], reverse(Σ_[2:end-1]),dims=1)
            ΣΣ_ = cat(reverse(Σ_[2:end-1]),Σ_[2:end-1],dims=1)
            #ΦΦ_ = cat(Φ_[2:end-1],-reverse(Φ_[2:end-1]),dims=1)
            ΦΦ_ = cat(-reverse(Φ_[2:end-1]),Φ_[2:end-1],dims=1)
            #Plots.plot!(xx_./xv,ΣΣ_./maximum(ΣΣ_),label="t: $tx",legend=false,xlim=(0,2*r.X●[end]/xv),ylim=(0,1.2))#*maximum(ΣΣ_)))#:bottomright)#,title="Σ")#;legend=false)
            #Plots.plot!(xx_./xv,ΣΣ_./maximum(ΣΣ_),label="t: $tx",legend=false,xlim=(-r.X●[end]/xv,r.X●[end]/xv),ylim=(0,1.2))#*maximum(ΣΣ_)))#:bottomright)#,title="Σ")#;legend=false)        display(fig)
            #Plots.plot!(xx_./xv,ΣΣ_/Σ_[end],xlim=(-3,3),ylim=(0,1.1*maximum(ΣΣ_./Σ_[end])), label=L"\log_{10} X = "*"$log10X",legend=:topright,linestyle=style,linewidth=1.5)#,title="Σ")#;legend=false)        display(fig)
            #Plots.plot!(xx_./xv,ΣΣ_/Σ_[end],ylim=(0,1.1*maximum(ΣΣ_./Σ_[end])), label=L"\log_{10} X = "*"$log10X",legend=:topright)#,linestyle=style,linewidth=1.5)#,title="Σ")#;legend=false)        display(fig)
            #Plots.plot!(xx_./xv,ΣΣ_,ylim=(0,1.1*maximum(ΣΣ_)), label=L"\log_{10} X = "*"$log10X",legend=:topright)#,linestyle=style,linewidth=1.5)#,title="Σ")#;legend=false)        display(fig)
            if PLOT
                Plots.plot!(xx_./xv,ΣΣ_)
                Plots.plot!(xx_./xv,ΦΦ_)
                print("time: $tx\n")
                #Plots.plot!(xx_./xv,ΦΦ_./maximum(ΦΦ_),label=false,linestyle=:dash)#,label="t: $t",legend=:topleft)#,title="J")#;legend=false)
                display(fig)
            end
        #end
        #    for j in 1:800
        #        Godunov.step!(r)#,dt)
        #    end
        Godunov.tstep!(r,dt)
    end

    # test
    #(1) Φ(xc) = -√q Σ(xc)?
    x_ = ((r.X●[1:end-1]+r.X●[2:end])/2.0)
    dx1 = x_[2]-x_[1]; dxe = x_[end]-x_[end-1]
    x_ = vcat(x_[1]-dx1,x_,x_[end]+dxe)
    Σinterp = linear_interpolation(x_,r.Σ_)
    Φinterp = linear_interpolation(x_,r.Φ_)
    Σc = Σinterp(xh)
    Φc = Φinterp(xh)
    sq = √(3ν/s)
    relerr = (Σc+sq*Φc)/Σc
    print("Σc: $Σc\n")
    print("Φc: $Φc\n")
    print("√q: $sq\n")
    print("(Σc + √q Φc)/Σc: $relerr\n")

    tf=r.t
    print("Done at t= $tf\n")
end