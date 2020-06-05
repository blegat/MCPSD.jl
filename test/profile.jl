using MCPSD
function g(L, n)
    for i in 1:n
        MCPSD.mc_psd(L, digits=12, iteration_limit=20, silent=true)
    end
end
using ProfileView
function f(sparse_L, T=Float64)
    A = T[
        0 1 0 0 1
        1 0 1 0 1
        0 1 0 1 0
        0 0 1 0 1
        1 1 0 1 0
    ]
    if sparse_L
        A = sparse(A)
    end
    L = Symmetric(-A)
    @profview g(L, 1000)
end
