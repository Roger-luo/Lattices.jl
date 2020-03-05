export Chain

struct Chain{BC, L}  <: BoundedLattice{BC, 1} end
Chain{BC}(n::Int) where BC = Chain{BC, n}()
Chain(n::Int) = Chain{Periodic}(n)

Base.length(::Chain{BC, L}) where {BC, L} = L
Base.size(::Chain{BC, L}) where {BC, L} = (L, )
Base.show(io::IO, x::Chain{BC, N}) where {BC, N} = print(io, "Chain{$BC}($N)")

to_coordinate(x::Chain, k) = Coordinate(k)
to_siteid(l::Chain, x::Coordinate{1}) = x.xs[1]

Base.size(::SiteIt{<:Chain{BC, N}}) where {BC, N} = (N, )

Base.length(::EdgeIt{<:Chain{Periodic, L}}) where L = L

function Base.iterate(::EdgeIt{<:Chain{Periodic, L}, K}, st=1) where {L, K}
    st > L && return nothing
    return (st, mod1(st + K, L)), st + 1
end
