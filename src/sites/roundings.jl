function Base.round(::Type{Bit{T}}, x::T) where {T <: Integer}
    x == zero(T) ? Bit(zero(T)) : Bit(one(T))
end

function Base.round(::Type{Bit{T}}, x::T; threshold::T=T(1e-3)) where {T <: AbstractFloat}
    x < threshold ? Bit(zero(T)) : Bit(one(T))
end

function Base.round(::Type{Spin{T}}, x::T) where T
    x > 0 ? Spin(one(T)) : Spin(-one(T))
end

function Base.round(::Type{Half{T}}, x::T) where {T <: AbstractFloat}
    x > 0 ? Half(T(0.5)) : Half(-T(0.5))
end
