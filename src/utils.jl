# Similar to `LinearAlgebra.isposdef`
function is_add_mul_psd(X::Symmetric, α, dX::Union{Diagonal, Symmetric})
    # With `check = true`, it throws an error,
    # with `check = false`, it does not and we can check with `issuccess`.
    return issuccess(cholesky!(sym_plus(X, _prod(α, dX))::typeof(X), check = false))
end

_prod(A, B) = A * B
@static if VERSION < v"1.4"
    # Needed for Julia v1.0, otherwise, it converts to Matrix
    _prod(A::Number, B::Symmetric) = Symmetric(A * parent(B), LinearAlgebra.sym_uplo(B.uplo))
end

# FIXME in Base, Symmetric{T, <:Diagonal} + Symmetric{T, <:SparseMatrixCSC} gives Symmetric{T, <:Array}
function sym_minus(A::Diagonal, B::Symmetric)
    Symmetric(A - parent(B), LinearAlgebra.sym_uplo(B.uplo))
end
@static if VERSION < v"1.4"
    sym_plus(A::Symmetric, B::Symmetric) = Symmetric(parent(A) + parent(B), LinearAlgebra.sym_uplo(B.uplo))
else
    sym_plus(A::Symmetric, B::Symmetric) = A + B
end
function sym_plus(A::Symmetric, B::Diagonal)
    Symmetric(parent(A) + B, LinearAlgebra.sym_uplo(A.uplo))
end

function _inv(A::Symmetric{<:Any, <:SparseMatrixCSC})
    D = Matrix(parent(A))
    @assert issymmetric(D)
    return inv(Symmetric(D))
end
_inv(A::Symmetric) = inv(A)
