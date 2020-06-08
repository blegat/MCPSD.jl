using MCPSD
using LinearAlgebra, SparseArrays, Test

import MathOptInterface
const MOI = MathOptInterface

function vect_upper(A)
    n = LinearAlgebra.checksquare(A)
    return [A[i, j] for j in 2:n for i in 1:(j - 1)]
end

function direct_test(L::Matrix{T}, expected_X, expected_y, expected_obj, tol) where T
    @testset "$(sparse_L ? "sparse" : "dense")" for sparse_L in [true, false]
        if sparse_L
            sL = Symmetric(sparse(parent(L)))
        else
            sL = Symmetric(L)
        end
        X, y, results = MCPSD.mc_psd(sL, digits=12, iteration_limit=20, silent=true)
        @test X isa Symmetric{T, Matrix{T}}
        @test X ≈ expected_X atol=tol rtol=tol
        @test y isa Vector{T}
        @test y ≈ expected_y atol=tol rtol=tol
        @test results[MCPSD.NumberOfIterations()] == 19
        @test results[MOI.SolveTime()] isa Float64
        @show results[MOI.SolveTime()]
    end
end

function wikipedia_example(T::Type)
    A = T[
        0 1 0 0 1
        1 0 1 0 1
        0 1 0 1 0
        0 0 1 0 1
        1 1 0 1 0
    ]
    L = -A
    expected_X = [
        1.0       -0.366823   0.12487    0.12487   -0.366823
       -0.366823   1.0       -0.968815   0.877204  -0.730881
        0.12487   -0.968815   1.0       -0.968815   0.877204
        0.12487    0.877204  -0.968815   1.0       -0.968815
       -0.366823  -0.730881   0.877204  -0.968815   1.0
    ]
    expected_y = [
        0.733646,
        2.066519,
        1.937629,
        1.937629,
        2.066519
    ]
    expected_obj = 8.741944118
    tol = 1e-5
    return -A, expected_X, expected_y, expected_obj, tol
end

@testset "Wikipedia example with $T" for T in [Float64, BigFloat]
    direct_test(wikipedia_example(T)...)
end

include("MOI_wrapper.jl")
