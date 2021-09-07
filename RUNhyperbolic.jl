include("godunov.jl")
using .Godunov
r=Godunov.ring()
Godunov.getAlphas!(r)
for i in 1:100
    Godunov.step!(r)
end
