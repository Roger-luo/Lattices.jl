export HilbertSpace

"""
    HilbertSpace{L, S, M}

`M` dimension Hilbert space of Site type `L` in with shape `S`.
"""
struct HilbertSpace{L, Shape, N, NSite}
    cache::Array{L, N}

    function HilbertSpace{L}(dims::Int...) where L
        Shape = Tuple{dims...}
        N = length(dims)
        NSite = prod(dims)
        new{L, Shape, N, NSite}(downs(L, dims))
    end
end

StaticArrays.Size(::Type{HilbertSpace{L, Shape}}) where {L, Shape} = Size(Shape)
Base.@pure Base.size(::Type{<:HilbertSpace{L, S}}) where {L, S} = tuple(S.parameters...)

function carrybit!(a::Array{L}) where L
    @inbounds for i in eachindex(a)
        if a[i] == up(L)
            a[i] = down(L)
        else
            a[i] = up(L)
            break
        end
    end
    a
end

Base.iterate(it::HilbertSpace{L}) where L = (fill!(it.cache, down(L)); iterate(it, 1))

function Base.iterate(it::HilbertSpace{L, S, N, NSite}, state) where {L, S, N, NSite}
    if_stop(it, state) && return nothing
    state == 1 && return (it.cache, state + 1)
    carrybit!(it.cache), state + 1
end

@generated function if_stop(it::HilbertSpace{L, S, N, NSite}, state) where {L, S, N, NSite}
    if sizeof(Int) * 8 > NSite
        :(state == $((1 << NSite) + 1))
    else
        quote
            flag = true
            for each in it.cache
                if each != up($L)
                    flag = false
                    break
                end
            end
            flag
        end
    end
end

Base.eltype(::HilbertSpace{L, S, N}) where {L, S, N} = Array{L, N}
Base.length(::HilbertSpace{L, S, N, NSite}) where {L, S, N, NSite} = 1 << NSite

function Base.collect(space::HilbertSpace)
    r = Vector{eltype(space)}(undef, length(space))
    @inbounds for (i, each) in enumerate(space)
        r[i] = copy(each)
    end
    r
end
