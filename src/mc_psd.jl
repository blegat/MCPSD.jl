function print_current(iter, cholcnt, αp, αd, δ, ψ, ϕ, start, current)
    time_sec = (current - start) / 1e9
    println(" time [s] | iter | #chol |   αp  |   αd  |  log(gap)  |   lower    |   upper")
    @printf(" %8.3f | %4.0d | %5.0d | %5.3f | %5.3f | %10.5f | %10.3f | %10.3f\n", time_sec, iter, cholcnt, αp, αd, log10(δ), ψ, ϕ)
end

function line_search(X::Symmetric{T}, ΔX::Union{Diagonal{T}, Symmetric{T}}) where T
    α = one(T)
    num_chol = 0
    while !is_add_mul_psd(X, α, ΔX)
        num_chol += 1
        α *= 0.8
    end
    # stay away from boundary
    if !isone(α)
        α *= 0.95
    end
    return α, num_chol + 1
end

"""
    mc_psd(L::Symmetric;
           digits=DEFAULT_VALUES[Digits()],
           silent::Bool=DEFAULT_VALUES[MOI.Silent()],
           time_limit_sec=DEFAULT_VALUES[MOI.TimeLimitSec()],
           iteration_limit=DEFAULT_VALUES[IterationLimit()]) where T

Solves the Max-Cut SDP problem:
```
max tr(L * X)     min sum(y)
    diag(X) = 1      Diag(y) ⪰ I
          X ⪰ 0
```
and returns the optimal values of `X`, `y` found, the number of iterations
and the time taken to solve. If `silent` is `true`, no output is printed.

Inspired by `mc_psd.m` written by Franz Rendl in February 1999.
"""
function mc_psd(L::Symmetric{T};
                digits=DEFAULT_VALUES[Digits()],
                silent::Bool=DEFAULT_VALUES[MOI.Silent()],
                time_limit_sec=DEFAULT_VALUES[MOI.TimeLimitSec()],
                iteration_limit=DEFAULT_VALUES[IterationLimit()]) where T

    start = time_ns()
    time_limit_ns = time_limit_sec * 1e9

    # initialize data
    n = LinearAlgebra.checksquare(L)             # problem size
    b = ones(T, n)
    X = Symmetric(Matrix(Diagonal(b)))            # (6.3)
    y = reshape(sum(abs, L, dims=1), n) .+ one(T) # (6.4) but not quite
    Z = sym_minus(Diagonal(y), L)::typeof(L)      # (6.5)
    ϕ = b ⋅ y
    ψ = L ⋅ X
    δ = ϕ - ψ

    μ = (Z ⋅ X) / (4 * n)
    αp = 1
    αd = 1
    iter = 0
    cholcnt = 0

    current = time_ns()
    silent || print_current(iter, cholcnt, αp, αd, δ, ψ, ϕ, start, current)

    gap_condition() = δ > max(abs(ϕ), one(T)) * (10one(T))^(-digits)

    while current - start < time_limit_ns && # while time limit not exceeded,
        iter < iteration_limit &&            # iteration limit not exceeded
        gap_condition()                      # and duality gap too large

        Zi::Symmetric{T} = _inv(Z)
        iter = iter + 1
        Δy =  (Zi .* X) \ (μ * diag(Zi) - b) # (6.8)
        ΔZ = Diagonal(Δy)                    # (6.9)
        ΔX = μ * Zi - X - Zi * ΔZ * X        # (6.10)
        ΔX = Symmetric((ΔX + ΔX') / 2)       # (6.11)

        # find steplengths αp and αd
        αp, _num_chol = line_search(X, ΔX)
        cholcnt += _num_chol
        αd, _num_chol = line_search(Z, ΔZ)
        cholcnt += _num_chol

        # update
        X = sym_plus(X, _prod(αp, ΔX))::Symmetric{T, Matrix{T}}
        y = y + αd * Δy
        Z = sym_plus(Z, αd * ΔZ)::typeof(L)
        μ = (X ⋅ Z) / (2 * n)
        if min(αp, αd) < 0.5
            μ = μ * 1.5
        end
        if αp + αd > 1.6
            μ = μ * 0.75
        end
        if αp + αd > 1.9
            # reduce μ, if step size good
            μ = μ / (1 + 0.1 * iter)
        end
        ϕ = b ⋅ y
        ψ = L ⋅ X
        δ = ϕ - ψ

        # display current iteration
        current = time_ns()
        silent || print_current(iter, cholcnt, αp, αd, δ, ψ, ϕ, start, current)
    end            # end of main loop

    results = Dict{MOI.AbstractModelAttribute, Any}(
        MOI.SolveTime() => (current - start) / 1e9,
        NumberOfIterations() => iter,
        NumberOfCholeskyFactorizations() => cholcnt,
    )

    if !gap_condition()
        results[MOI.TerminationStatus()] = MOI.OPTIMAL
        results[MOI.PrimalStatus()] = MOI.FEASIBLE_POINT
        results[MOI.DualStatus()] = MOI.FEASIBLE_POINT
    elseif iter >= iteration_limit
        results[MOI.TerminationStatus()] = MOI.ITERATION_LIMIT
        results[MOI.PrimalStatus()] = MOI.UNKNOWN_RESULT_STATUS
        results[MOI.DualStatus()] = MOI.UNKNOWN_RESULT_STATUS
    else
        results[MOI.TerminationStatus()] = MOI.TIME_LIMIT
        results[MOI.PrimalStatus()] = MOI.UNKNOWN_RESULT_STATUS
        results[MOI.DualStatus()] = MOI.UNKNOWN_RESULT_STATUS
    end

    return X, y, results
end
