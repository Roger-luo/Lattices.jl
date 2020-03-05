export Boundary, Open, Periodic, sites, edges, faces, Coordinate

abstract type Boundary end
struct Open <: Boundary end
struct Periodic <: Boundary end

abstract type AbstractLattice{N} end
abstract type BoundedLattice{B <: Boundary, N} <: AbstractLattice{N} end

Base.ndims(::AbstractLattice{N}) where N = N

struct SiteIt{L <: AbstractLattice}
    lattice::L
end

Base.length(it::SiteIt) = length(it.lattice)
Base.eltype(::SiteIt) = Int

function Base.iterate(it::SiteIt, st=1)
    st > length(it) && return nothing
    st, st + 1
end

struct EdgeIt{L <: AbstractLattice, K}
    lattice::L
end

Base.eltype(::EdgeIt) = Tuple{Int, Int}

struct FaceIt{L <: AbstractLattice{2}, K}
    lattice::L
end

Base.eltype(::FaceIt) = NTuple{4, Int}

sites(x) = SiteIt(x)
edges(x; distance=1) = EdgeIt{typeof(x), distance}(x)
faces(x; distance=1) = FaceIt{typeof(x), distance}(x)

Base.show(io::IO, x::SiteIt) = print(io, "sites(", x.lattice, ")")
Base.show(io::IO, x::EdgeIt{<:Any, K}) where K = print(io, "edges(", x.lattice, "; distance=", K, ")")
Base.show(io::IO, x::FaceIt{<:AbstractLattice{2}, K}) where K = print(io, "faces(", x.lattice, "; distance=", K, ")")


struct Coordinate{N}
    xs::NTuple{N, Int}
end

Coordinate(xs...) = Coordinate(xs)

Base.length(::Coordinate{N}) where N = N
Base.size(::Coordinate{N}) where N = (N, )
Base.iterate(coo::Coordinate, st=1) = iterate(coo.xs, st)
Base.getindex(coo::Coordinate, idx::Int) = coo.xs[idx]
Base.eltype(::Coordinate) = Int
Base.show(io::IO, x::Coordinate) = print(io, "Coordinate", x.xs)
Base.show(io::IO, x::Coordinate{1}) = print(io, "Coordinate(", x.xs[1], ")")
