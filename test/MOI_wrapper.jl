using Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

import MCPSD

@testset "SolverName" begin
    @test MOI.get(MCPSD.Optimizer(), MOI.SolverName()) == "MCPSD"
end

@testset "supports_incremental_interface" begin
    @test MOI.supports_incremental_interface(MCPSD.Optimizer())
end

function moi_test(optimizer, L::Matrix{T}, expected_X, expected_y, expected_obj, tol) where T
    @test MOI.is_empty(optimizer)
    @test MOI.get(optimizer, MOI.SolverName()) == "MCPSD"
    MOI.supports(optimizer, MCPSD.Digits())
    MOI.set(optimizer, MCPSD.Digits(), 12)
    @test MOI.get(optimizer, MCPSD.Digits()) == 12
    MOI.supports(optimizer, MCPSD.IterationLimit())
    MOI.set(optimizer, MCPSD.IterationLimit(), 20)
    @test MOI.get(optimizer, MCPSD.IterationLimit()) == 20
    MOI.supports(optimizer, MOI.Silent())
    MOI.set(optimizer, MOI.Silent(), true)
    @test MOI.get(optimizer, MOI.Silent())
    @test MOI.is_empty(optimizer)
    @test isempty(MOI.get(optimizer, MOI.ListOfConstraintTypesPresent()))
    F = MOI.VectorOfVariables
    S = MCPSD.Elliptope
    @test MOI.get(optimizer, MOI.NumberOfConstraints{F, S}()) == 0
    @test isempty(MOI.get(optimizer, MOI.ListOfConstraintIndices{F, S}()))
    x, cx = MOI.add_constrained_variables(optimizer, MCPSD.Elliptope(size(L, 1)))
    @test !MOI.is_empty(optimizer)
    @test cx == MOI.ConstraintIndex{F, S}(1)
    @test MOI.get(optimizer, MOI.ListOfConstraintTypesPresent()) == [(F, S)]
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

    MOI.get(optimizer, MOI.ResultCount()) == 0
    MOI.optimize!(optimizer)
    MOI.get(optimizer, MOI.ResultCount()) == 1
    expected_x = vect_upper(expected_X)
    @test MOI.get(optimizer, MOI.TerminationStatus()) == MOI.OPTIMAL
    @test MOI.get(optimizer, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    @test MOI.get(optimizer, MOI.DualStatus()) == MOI.FEASIBLE_POINT
    @test MOI.get(optimizer, MOI.ObjectiveValue()) ≈ expected_obj atol=tol rtol=tol
    @test MOI.get(optimizer, MOI.DualObjectiveValue()) ≈ expected_obj atol=tol rtol=tol
    @test MOI.get(optimizer, MOI.VariablePrimal(), x) ≈ expected_x atol=tol rtol=tol
    @test MOI.get(optimizer, MOI.ConstraintPrimal(), cx) ≈ expected_x atol=tol rtol=tol
    @test MOI.get(optimizer, MOI.ConstraintDual(), cx) ≈ expected_y atol=tol rtol=tol
    @test MOI.get(optimizer, MCPSD.NumberOfIterations()) == 19
    @test MOI.get(optimizer, MOI.SolveTimeSec()) isa Float64
    @show MOI.get(optimizer, MOI.SolveTimeSec())
end

@testset "[MOI_wrapper] Wikipedia example with $T" for T in [Float64, BigFloat]
    optimizer = MCPSD.Optimizer{T}()
    moi_test(optimizer, wikipedia_example(T)...)
    cached = MOIU.CachingOptimizer(MCPSD.Optimizer{T}(), MOIU.AUTOMATIC)
    MOI.empty!(optimizer)
    MOIU.reset_optimizer(cached, optimizer)
    moi_test(cached, wikipedia_example(T)...)
end
