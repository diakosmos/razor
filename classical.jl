module Classical
using LinearAlgebra
#using StaticArrays

"""
visco()

Implements a power-law dependence for the viscosity ν on the surface density Σ,

    ν = ν₀(Σ/Σ₀)ᵝ

following Grätz, Seiß & Spahn (2018).
"""
function visco(Σrel,β)
    THRESH = 0.01 # Thresholding allows us to use β<0, which would otherwise lead to infinities.
    threshold(Σ) = Σ > THRESH ? Σ : THRESH
    Σr = threshold.(Σrel)
    return Diagonal(Σr.^β)
end

"""
Potentially confusing terminology: a planetary particulate disk is called a "ring",
not a disk.
"""
mutable struct ring
#    α::Float64
#    Δx::Float64
#    Δt::Float64
#    xi, xf ::Float64
    N::Integer
    β::Float64
    ΔtMAX::Float64
    # Defined on nodes:
    X_::Array{Float64,1}
    ΔX_::Array{Float64,1}
    U_::Array{Float64,1}
    # Defined on cells:
    x_::Array{Float64,1}
    σ_::Array{Float64,1}
    Gᵐ::Any # Tridiagonal
#    μ_::Array{Float64,1}
#    νᵐ::Any
end
function ring(N=5;β=0,xf=3.0,xi=0.01)#10*xf/N)
    xtemp_ = collect(range(xi,xf,length=N+1))
    X_ = xtemp_[2:end]
    U_ = 1.0 ./ X_.^4
    ΔX_ = diff(xtemp_)
    x_ = 0.5*(xtemp_[1:end-1]+xtemp_[2:end]) # midpoints
    σ_ = ones(N)
    νᵐ = visco(σ_,β)
    Mᵐ = Diagonal(ΔX_)
    Mⁱᵐ= Diagonal(1.0./ΔX_)
    Kdiag = [1.0 ./ diff(x_)..., 0.0]
    Ksuper = -Kdiag[1:end-1]
    Ksub = zeros(N-1)
    Kᵐ = Tridiagonal(Ksub,Kdiag,Ksuper)
    Dᵐ = Tridiagonal(ones(N-1),-ones(N),zeros(N-1))
    Gᵐ = 3 * Dᵐ * Kᵐ * νᵐ
    μ_ = Mᵐ * σ_ # mass (per azimuthal length) in annulus
    ΔtVN = 0.5 * minimum(diff(X_).^2 ./ diag(νᵐ)[1:end-1] )
    ΔtADV = minimum(ΔX_ ./ U_)
    ΔtMAX = min(ΔtVN,ΔtADV)
    return ring(N,β,ΔtMAX,X_,ΔX_,U_,x_,σ_,Gᵐ)#,μ_)
end


function adjust!(r::ring)
    νᵐ = visco(r.σ_,r.β)
    Kdiag = [1.0 ./ diff(r.x_)..., 0.0]
    Ksuper = -Kdiag[1:end-1]
    Ksub = zeros(r.N-1)
    Kᵐ = Tridiagonal(Ksub,Kdiag,Ksuper)
    Dᵐ = Tridiagonal(ones(r.N-1),-ones(r.N),zeros(r.N-1))
    Gᵐ = 3 * Dᵐ * Kᵐ * νᵐ
    r.Gᵐ = Gᵐ
end


function advect!(r::ring;Δt=r.ΔtMAX)
    σ = r.σ_; U = r.U_;; ΔX = r.ΔX_
    σ[2:end]   += Δt .* (σ[1:end-1] .* U[1:end-1] ./ ΔX[2:end] )
    σ[1:end] -= Δt .* (σ[1:end] .* U[1:end] ./ ΔX[1:end] )
end

function diffuse!(r::ring;Δt=r.ΔtMAX,α=0.5)
    Gᵐ = r.Gᵐ
    N = r.N
    Aᵐ = I(N) -    α  * Δt * Gᵐ
    Bᵐ = I(N) + (1-α) * Δt * Gᵐ
    # Solve:
    r.σ_ = Aᵐ \ (Bᵐ * r.σ_)
#    r.σ_[1]   = 0.0
    r.σ_[end] = 1.0
end

function iterate!(r::ring;Δt=r.ΔtMAX,α=0.5)
    adjust!(r)
    diffuse!(r;Δt=Δt,α=α)
    advect!(r;Δt=Δt)
end

end#module
