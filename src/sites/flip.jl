export flip, flip!, randflip!, swap!, randswap!, randswap

"""
    FlipStyle

Abstract type for flip styles.
"""
abstract type FlipStyle end

"""
    flip(site)

Defines how to flip this type of `site`.
"""
function flip end

flip(s::BT) where {BT <: BinarySite} = s == up(BT) ? down(BT) : up(BT)

"""
    flip!(S, I...)

Flips given configuration `S` at index `I...`.
"""
function flip! end

flip!(S::AbstractArray{BT, N}, I...) where {BT <: BinarySite, N} = flip!(S, I)
flip!(S::AbstractArray{BT, N}, I::Tuple) where {BT <: BinarySite, N} = (@inbounds S[I...] = flip(S[I...]); S)

"""
    randflip!(config) -> (proposed_index, config)
    randflip!(::FlipStyle, config) -> (proposed_index, config)

Flip the lattice configuration randomly using given flip style, default flip style is
[`UniformFlip`](@ref){1}, which choose **one** site in the configuration uniformly and flip it.

One should always be able to use `config[proposed_index]` to get the current value of
this lattice configuration. Whether one can change the site, depends on whether the
configuration is stored in a mutable type.
"""
function randflip! end

randflip!(S::AbstractArray) = randflip!(FlipStyle(S), S)

"""
    UniformFlip{N}

Choose `N` sites in the configuration uniformly and flip it
"""
struct UniformFlip{N} <: FlipStyle end
UniformFlip(N::Int) = UniformFlip{N}()
FlipStyle(S::AbstractArray{BT, N}) where {BT <: BinarySite, N} = UniformFlip(1) # use uniform(1) as default

function randflip!(::UniformFlip{1}, S::AbstractArray{BT}) where {BT <: BinarySite}
    proposed_index = rand(1:length(S))
    flip!(S, proposed_index)
    proposed_index, S
end

# TODO:
# 1. minimum sampling length of index can be calculated
#    use this to generate a shorter sequence
# 2. specialized method for static arrays (their length can be found in compile time)
function randflip!(::UniformFlip{N}, S::AbstractArray{BT}) where {BT <: BinarySite, N}
    proposed = zeros(Int, N)

    count = 0
    while count != N
        i = rand(1:length(S))
        i in proposed || (count += 1; proposed[count] = i)
    end

    @inbounds for i in 1:N
        flip!(S, proposed[i])
    end

    proposed, S
end

function swap!(S::AbstractArray, I1, I2)
    S[I2], S[I1] = S[I1], S[I2]
    return S
end

randswap!(S::AbstractArray) = randswap!(UniformFlip(1), S)

function randswap!(::UniformFlip{1}, S::AbstractArray)
    proposed = rand(1:length(S)), rand(1:length(S))
    swap!(S, proposed...)
    return proposed, S
end

randswap(S::AbstractArray) = randswap!(copy(S))
