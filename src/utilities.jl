@inline _float_type(::Type{T}) where {T<:Real} = float(T)

@inline _promote_float_type(x::Real) = _float_type(typeof(x))

@inline _promote_float_type(x::Real, y::Real) =
    promote_type(_promote_float_type(x), _promote_float_type(y))

@inline _promote_float_type(x::Real, y::Real, z::Real...) =
    promote_type(_promote_float_type(x, y), _promote_float_type(z...))

@inline _halfpi(::Type{T}) where {T<:AbstractFloat} = T(π / 2)
@inline _pi(::Type{T}) where {T<:AbstractFloat} = T(π)
