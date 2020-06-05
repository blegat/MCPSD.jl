struct Elliptope <: MOI.AbstractVectorSet
    side_dimension::Int
end
MOI.dimension(set::Elliptope) = MOI.dimension(MOI.PositiveSemidefiniteConeTriangle(set.side_dimension - 1))

struct Digits <: MOI.AbstractOptimizerAttribute end
struct IterationLimit <: MOI.AbstractOptimizerAttribute end

const SUPPORTED_SETTABLE_ATTRIBUTES = Union{
    Digits,
    MOI.Silent,
    MOI.TimeLimitSec,
    IterationLimit
}
const DEFAULT_VALUES = Dict{MOI.AbstractOptimizerAttribute, Any}(
    Digits() => 5.5,
    MOI.Silent() => false,
    MOI.TimeLimitSec() => 10.0,
    IterationLimit() => 100
)

struct NumberOfIterations <: MOI.AbstractModelAttribute
end
MOI.is_set_by_optimize(::NumberOfIterations) = true

struct NumberOfCholeskyFactorizations <: MOI.AbstractModelAttribute
end
MOI.is_set_by_optimize(::NumberOfCholeskyFactorizations) = true

const SUPPORTED_GETTABLE_ATTRIBUTES = Union{
    MOI.TerminationStatus,
    MOI.PrimalStatus,
    MOI.DualStatus,
    MOI.SolveTime,
    NumberOfIterations,
    NumberOfCholeskyFactorizations
}
