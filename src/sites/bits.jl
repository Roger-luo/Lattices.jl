struct FixedLengthBits{L, N}
    chunks::NTuple{N, UInt64}
end

FixedLengthBits(x::Int) = FixedLengthBits{64 - leading_zeros(x), 1}((UInt(x), ))

const _msk64 = ~UInt64(0)
@inline _div64(l) = l >> 6
@inline _mod64(l) = l & 63
@inline _msk_end(l::Integer) = _msk64 >>> _mod64(-l)
@inline _msk_end(B::FixedLengthBits{L}) where L = _msk_end(L)
num_bit_chunks(n::Int) = _div64(n+63)

function Base.getindex(x::FixedLengthBits{L}, i::Int) where L
    @boundscheck 0 < i ≤ L ? true : throw(BoundsError(x, i))
    Int((x.chunks[_div64(i) + 1] >> (_mod64(i) - 1)) & 0x01)
end

function Base.iterate(x::FixedLengthBits{L}, state=1) where L
    state == L + 1 && return nothing
    @inbounds(x[state]), state + 1
end

Base.length(::FixedLengthBits{L}) where L = L
Base.eltype(::FixedLengthBits) = Int

function Base.show(io::IO, x::FixedLengthBits{L}) where L
    @inbounds for i = L:-1:1
        print(io, x[i])
    end
end
