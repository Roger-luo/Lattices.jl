module Lattices

include("base.jl")

# lib
include("lib/chain.jl")
include("lib/square.jl")

# sites
using Random, StaticArrays
include("sites/spins.jl")
include("sites/space.jl")
include("sites/flip.jl")

ups(::Type{ST}, ltc::AbstractLattice) where ST = ups(ST, size(ltc))
downs(::Type{ST}, ltc::AbstractLattice) where ST = downs(ST, size(ltc))
Random.rand(::Type{ST}, ltc::AbstractLattice) where ST = rand(ST, size(ltc))

ups(ltc::AbstractLattice) = ups(Bit{Int}, ltc)
downs(ltc::AbstractLattice) = downs(Bit{Int}, ltc)
Random.rand(ltc::AbstractLattice) = rand(Bit{Int}, ltc)

HilbertSpace{L}(ltc::AbstractLattice) where L = HilbertSpace{L}(size(ltc)...)
HilbertSpace(ltc::AbstractLattice) = HilbertSpace{Bit{Int}}(ltc)

export SU
struct SU{N} end
SU(n::Int) = SU{n}()
Base.show(io::IO, x::SU{N}) where N = print(io, "SU(", N, ")")

Random.rand(s::SU{2}, ltc::AbstractLattice) = rand(Bit{Int}, s, ltc)
function Random.rand(::Type{T}, ::SU{2}, ltc::AbstractLattice) where T
    L = length(ltc)
    indices = randperm(L)[1:L÷2]
    S = downs(T, ltc)
    for each in indices
        S[each] = up(T)
    end
    return S
end

end # module
