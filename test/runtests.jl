using MCPSD
using LinearAlgebra, SparseArrays, Test

import MathOptInterface
const MOI = MathOptInterface

function vect_upper(A)
    n = LinearAlgebra.checksquare(A)
    return [A[i, j] for j in 2:n for i in 1:(j - 1)]
end

function direct_test(L::Matrix{T}, expected_X, expected_y, tol) where T
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

function moi_test(L::Matrix{T}, expected_X, expected_y, tol) where T
    optimizer = MCPSD.Optimizer{T}()
    MOI.set(optimizer, MCPSD.Digits(), 12)
    @test MOI.get(optimizer, MCPSD.Digits()) == 12
    MOI.set(optimizer, MCPSD.IterationLimit(), 20)
    @test MOI.get(optimizer, MCPSD.IterationLimit()) == 20
    MOI.set(optimizer, MOI.Silent(), true)
    @test MOI.get(optimizer, MOI.Silent())
    @test isempty(MOI.get(optimizer, MOI.ListOfConstraints()))
    F = MOI.VectorOfVariables
    S = MCPSD.Elliptope
    @test MOI.get(optimizer, MOI.NumberOfConstraints{F, S}()) == 0
    @test isempty(MOI.get(optimizer, MOI.ListOfConstraintIndices{F, S}()))
    x, cx = MOI.add_constrained_variables(optimizer, MCPSD.Elliptope(size(L, 1)))
    @test cx == MOI.ConstraintIndex{F, S}(1)
    @test MOI.get(optimizer, MOI.ListOfConstraints()) == [(F, S)]
    @test MOI.get(optimizer, MOI.NumberOfConstraints{F, S}()) == 1
    @test MOI.get(optimizer, MOI.ListOfConstraintIndices{F, S}()) == [cx]
    @test all(x -> MOI.is_valid(optimizer, x), x)
    @test MOI.is_valid(optimizer, cx)
    l = vect_upper(L)
    func = MOI.ScalarAffineFunction{T}(MOI.ScalarAffineTerm{T}[MOI.ScalarAffineTerm(2l[i], x[i]) for i in eachindex(x)], zero(T))
    @test MOI.get(optimizer, MOI.ObjectiveSense()) == MOI.FEASIBILITY_SENSE
    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    @test MOI.get(optimizer, MOI.ObjectiveSense()) == MOI.MAX_SENSE
    MOI.set(optimizer, MOI.ObjectiveFunction{typeof(func)}(), func)
    @test MOI.get(optimizer, MOI.ObjectiveFunction{typeof(func)}()) ≈ func
    MOI.optimize!(optimizer)
    expected_x = vect_upper(expected_X)
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test MOI.get(optimizer, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    @test MOI.get(optimizer, MOI.DualStatus()) == MOI.FEASIBLE_POINT
    @test MOI.get(optimizer, MOI.VariablePrimal(), x) ≈ expected_x atol=tol rtol=tol
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), cx) ≈ expected_x atol=tol rtol=tol
    @test MOI.get(optimizer, MOI.ConstraintDual(), cx) ≈ expected_y atol=tol rtol=tol
    @test MOI.get(optimizer, MCPSD.NumberOfIterations()) == 19
    @test MOI.get(optimizer, MOI.SolveTime()) isa Float64
    @show MOI.get(optimizer, MOI.SolveTime())
end

@testset "Wikipedia example with $T" for T in [Float64, BigFloat]
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
    tol = 1e-5
    expected_y = [
        0.733646,
        2.066519,
        1.937629,
        1.937629,
        2.066519
    ]
    direct_test(L, expected_X, expected_y, tol)
    moi_test(L, expected_X, expected_y, tol)
end
