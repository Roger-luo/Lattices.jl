export up, down, Bit, Spin, Half, Clock, Potts, Continuous

abstract type AbstractSite{T} end

Base.show(io::IO, s::AbstractSite) = show(io, value(s))

Base.:(==)(lhs::AbstractSite, rhs::Number) = value(lhs) == rhs
Base.:(==)(lhs::Number, rhs::AbstractSite) = lhs == value(rhs)

abstract type IntegerSite{T} <: AbstractSite{T} end
abstract type BinarySite{T} <: IntegerSite{T} end
abstract type ContinuousSite{T} <: AbstractSite{T} end

"""
    Bit{T} <: BinarySite{T}
"""
struct Bit{T} <: BinarySite{T}
    value::T
end

struct Spin{T} <: BinarySite{T}
    value::T
end

struct Half{T} <: BinarySite{T}
    value::T
end

struct Clock{T, q} <: IntegerSite{T}
    value::T
end

struct Potts{T, q} <: IntegerSite{T}
    value::T
end

struct Continuous{T, d} <: ContinuousSite{T}
    value::T
end

"""
    up(site) -> site
    up(site_type) -> site

Up (highest value) tag for this label. e.g. `1` for `Bit`, `0.5` for `Half`.
"""
up(::ST) where {ST<: IntegerSite} = up(ST)
up(::Type{ST}) where {ST <: IntegerSite} = values(ST)[end]

"""
    down(site) -> site
    down(site_type) -> site

Down (lowest value) tag for this label. e.g. `0` for `Bit`, `-0.5` for `Half`.
"""
down(::ST) where {ST<: IntegerSite} = down(ST)
down(::Type{ST}) where {ST <: IntegerSite} = values(ST)[1]

"""
    values(site) -> site
    values(site_type) -> site

Returns a tuple of all possible values of the site type
"""
Base.values(::ST) where {ST <: AbstractSite} = values(ST)
Base.values(::Type{ST}) where {ST <: AbstractSite} = map(ST, _values(ST))
Base.values(::ST, i::Int) where {ST <: AbstractSite} = _values(ST)[i]

_values(::Type{Bit{T}}) where T = (zero(T), one(T))
_values(::Type{Spin{T}}) where T = (-one(T), one(T))
_values(::Type{Half{T}}) where T = (-0.5, 0.5)
_values(::Type{Clock{T, q}}) where {T, q} = Base.OneTo(q)
_values(::Type{Potts{T, q}}) where {T, q} = Tuple(-q:q)
# _values(::Type{Continuous{T}}) where {T} = # TODO

value(S::AbstractSite) = S.value

Base.length(::AbstractSite) = 1
Base.iterate(x::AbstractSite) = (x, nothing)
Base.iterate(x::AbstractSite, state) = nothing

Random.rand(rng::Random.AbstractRNG, sp::Random.SamplerType{Bit{T}}) where T = Bit{T}(rand(rng, Bool))
Random.rand(rng::Random.AbstractRNG, sp::Random.SamplerType{Spin{T}}) where T = Spin{T}(2 * rand(rng, Bool) - 1)
Random.rand(rng::Random.AbstractRNG, sp::Random.SamplerType{Half{T}}) where T = Half{T}(rand(rng, Bool) - 0.5)
Random.rand(rng::Random.AbstractRNG, sp::Random.SamplerType{Clock{T, q}}) where {T, q} = Clock{T, q}(rand(rng, 1:q))
Random.rand(rng::Random.AbstractRNG, sp::Random.SamplerType{Potts{T, q}}) where {T, q} = Potts{T, q}(rand(rng, -q:q))

Base.to_index(A::AbstractArray, i::Bit) = Int(value(i) + 1)
Base.to_index(A::AbstractArray, i::Spin) = Int(0.5 * (value(i) + 1) + 1)
Base.to_index(A::AbstractArray, i::Half) = Int(value(i) + 1.5)
Base.to_index(A::AbstractArray, i::Clock) = Int(value(i))
Base.to_index(A::AbstractArray, i::Potts) = Int(value(i) + q + 1)

Base.to_index(A::AbstractArray, I::AbstractArray{Bit{T}, N}) where {T, N} = Int(I) + 1

Base.getindex(t::Tuple, i::Bit) = getindex(t, value(i) + 1)
Base.getindex(t::Tuple, i::Spin) = getindex(t, Int(div(value(i) + 1, 2) + 1))

include("roundings.jl")
include("conversions.jl")
include("arraymath.jl")
