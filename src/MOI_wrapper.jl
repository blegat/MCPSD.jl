mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    sense::MOI.Union{Nothing, MOI.OptimizationSense}
    objective_function::Union{Nothing, MOI.ScalarAffineFunction{T}}
    set::Union{Nothing, Elliptope}
    X::Union{Nothing, Symmetric{T, Matrix{T}}}
    x::Union{Nothing, Vector{T}}
    y::Union{Nothing, Vector{T}}
    results::Union{Nothing, Dict{MOI.AbstractModelAttribute, Any}}
    options::Dict{MOI.AbstractOptimizerAttribute, Any}
    function Optimizer{T}() where T
        return new{T}(nothing, nothing, nothing, nothing, nothing, nothing, nothing,
                   Dict{MOI.AbstractOptimizerAttribute, Any}())
    end
end
Optimizer() = Optimizer{Float64}()

MOI.get(::Optimizer, ::MOI.SolverName) = "MCPSD"
MOI.supports(::Optimizer, ::SUPPORTED_SETTABLE_ATTRIBUTES) = true
function MOI.set(optimizer::Optimizer, attr::SUPPORTED_SETTABLE_ATTRIBUTES, value)
    optimizer.options[attr] = value
end
function MOI.get(optimizer::Optimizer, attr::SUPPORTED_SETTABLE_ATTRIBUTES)
    return get(optimizer.options, attr, DEFAULT_VALUES[attr])
end

MOI.get(optimizer::Optimizer, ::MOI.ResultCount) = optimizer.results === nothing ? 0 : 1
function MOI.get(optimizer::Optimizer, attr::SUPPORTED_GETTABLE_ATTRIBUTES)
    optimizer.results === nothing && throw(MOI.OptimizeNotCalled())
    return optimizer.results[attr]
end

MOI.get(optimizer::Optimizer, attr::MOI.ObjectiveValue) = MOI.Utilities.get_fallback(optimizer, attr)
MOI.get(optimizer::Optimizer, ::MOI.DualObjectiveValue) = sum(optimizer.y)

function MOI.get(optimizer::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex)
    MOI.throw_if_not_valid(optimizer, vi)
    optimizer.x === nothing && throw(MOI.OptimizeNotCalled())
    return optimizer.x[vi.value]
end
function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintPrimal,
                 ci::MOI.ConstraintIndex{MOI.VectorOfVariables, Elliptope})
    return MOI.Utilities.get_fallback(optimizer, attr, ci)
end
function MOI.get(optimizer::Optimizer, ::MOI.ConstraintDual,
                 ci::MOI.ConstraintIndex{MOI.VectorOfVariables, Elliptope})
    MOI.throw_if_not_valid(optimizer, ci)
    optimizer.y === nothing && throw(MOI.OptimizeNotCalled())
    return optimizer.y
end

function MOI.is_empty(optimizer::Optimizer)
    return optimizer.sense === nothing && optimizer.objective_function === nothing &&
        optimizer.set === nothing
end

function MOI.empty!(optimizer::Optimizer)
    optimizer.sense = nothing
    optimizer.objective_function = nothing
    optimizer.set = nothing
    optimizer.X = nothing
    optimizer.x = nothing
    optimizer.y = nothing
    optimizer.results = nothing
end

MOI.Utilities.supports_default_copy_to(::Optimizer, names::Bool) = !names

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kwargs...)
    return MOI.Utilities.automatic_copy_to(dest, src; kwargs...)
end

function MOI.get(optimizer::Optimizer, ::MOI.ListOfModelAttributesSet)
    attributes = MOI.AbstractModelAttribute[]
    if optimizer.sense !== nothing
        push!(optimizer, optimizer.sense)
    end
    if optimizer.objective_function !== nothing
        F = typeof(optimizer.objective_function)
        push!(optimizer, MOI.ObjectiveFunction{F}())
    end
    return attributes
end

function MOI.get(::Optimizer, ::MOI.ListOfVariableAttributesSet)
    return MOI.AbstractVariableAttribute[]
end

function MOI.get(::Optimizer, ::MOI.ListOfConstraintAttributesSet)
    return MOI.AbstractConstraintAttribute[]
end

function MOI.supports(
    optimizer::Optimizer{T},
    ::Union{MOI.ObjectiveSense,
            MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}}) where T
    return true
end

MOI.supports_add_constrained_variables(::Optimizer, ::Type{<:Elliptope}) = true
function MOI.add_constrained_variables(optimizer::Optimizer, set::Elliptope)
    optimizer.set = set
    return MOI.VariableIndex.(1:MOI.dimension(set)),
        MOI.ConstraintIndex{MOI.VectorOfVariables, Elliptope}(1)
end
function MOI.is_valid(optimizer::Optimizer, vi::MOI.VariableIndex)
    return optimizer.set !== nothing && 1 <= vi.value <= MOI.dimension(optimizer.set)
end
function MOI.is_valid(optimizer::Optimizer, ci::MOI.ConstraintIndex{MOI.VectorOfVariables, Elliptope})
    return optimizer.set !== nothing && isone(ci.value)
end
function MOI.get(optimizer::Optimizer, ::MOI.ConstraintFunction,
                 ci::MOI.ConstraintIndex{MOI.VectorOfVariables, Elliptope})
    MOI.throw_if_not_valid(optimizer, ci)
    return MOI.VectorOfVariables(MOI.get(optimizer, MOI.ListOfVariableIndices()))
end
function MOI.get(optimizer::Optimizer, ::MOI.NumberOfVariables)
    return optimizer.set === nothing ? 0 : MOI.dimension(optimizer.set)
end
function MOI.get(optimizer::Optimizer, ::MOI.ListOfVariableIndices)
    return MOI.VariableIndex.(1:MOI.get(optimizer, MOI.NumberOfVariables()))
end
MOI.get(optimizer::Optimizer, ::MOI.VariableName, ::MOI.VariableIndex) = ""
function MOI.get(optimizer::Optimizer, ::MOI.ListOfConstraints)
    constraint_types = Tuple{DataType, DataType}[]
    if optimizer.set !== nothing
        push!(constraint_types, (MOI.VectorOfVariables, Elliptope))
    end
    return constraint_types
end
function MOI.get(optimizer::Optimizer, ::MOI.NumberOfConstraints{MOI.VectorOfVariables, Elliptope})
    return optimizer.set === nothing ? 0 : 1
end
function MOI.get(optimizer::Optimizer, ::MOI.ListOfConstraintIndices{MOI.VectorOfVariables, Elliptope})
    indices = MOI.ConstraintIndex{MOI.VectorOfVariables, Elliptope}[]
    if optimizer.set !== nothing
        push!(indices, MOI.ConstraintIndex{MOI.VectorOfVariables, Elliptope}(1))
    end
    return indices
end

function MOI.set(optimizer::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    if sense == MOI.FEASIBILITY_SENSE
        optimizer.objective_function = nothing
    end
    optimizer.sense = sense
end
function MOI.get(optimizer::Optimizer, ::MOI.ObjectiveSense)
    return optimizer.sense === nothing ? MOI.FEASIBILITY_SENSE : optimizer.sense
end
function MOI.set(optimizer::Optimizer{T}, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}},
                 func::MOI.ScalarAffineFunction{T}) where T
    optimizer.objective_function = func
end
function MOI.get(optimizer::Optimizer{T}, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}) where T
    return optimizer.objective_function
end
function MOI.optimize!(optimizer::Optimizer{T}) where T
    !iszero(optimizer.objective_function.constant) && error("Nonzero objective constant not supported")
    mapping = [(i, j) for j in 2:optimizer.set.side_dimension for i in 1:(j - 1)]
    I = Int[]
    J = Int[]
    V = T[]
    for term in optimizer.objective_function.terms
        i, j = mapping[term.variable_index.value]
        push!(I, i)
        push!(I, j)
        push!(J, j)
        push!(J, i)
        coef = term.coefficient
        coef = (optimizer.sense == MOI.MIN_SENSE ? -coef : coef) / 2
        push!(V, coef)
        push!(V, coef)
    end
    L = Symmetric(sparse(I, J, V))
    optimizer.X, optimizer.y, optimizer.results = mc_psd(L,
       digits = MOI.get(optimizer, Digits()),
       silent = MOI.get(optimizer, MOI.Silent()),
       time_limit_sec = MOI.get(optimizer, MOI.TimeLimitSec()),
       iteration_limit = MOI.get(optimizer, IterationLimit()))
    optimizer.x = [optimizer.X[i, j] for j in 2:optimizer.set.side_dimension for i in 1:(j - 1)]
    return
end
