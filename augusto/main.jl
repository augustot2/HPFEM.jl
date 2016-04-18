include("../src/HPFEM.jl")
include("lagrange_plot.jl")
using PyPlot
using HPFEM
using Interact
using GtkInteract


fun(x) = (1 + 4*pi^2)*sin(2*pi*x)
resp(x) = sin(2*pi*x)

f = figure()
@manipulate for Q= 1:1:100,M=1:1:100, nel=1:1:100 ; withfig(f) do
        if Q < M
            Q=M
        end
        lagrange_oed_plot(M,Q,nel,fun,resp)
    end
end
