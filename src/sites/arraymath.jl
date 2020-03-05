export ups, downs

ups(::Type{ST}, dims::Int...) where ST = ups(ST, dims)
downs(::Type{ST}, dims::Int...) where ST = downs(ST, dims)

ups(::Type{ST}, dims::Dims) where ST = fill(up(ST), dims)
downs(::Type{ST}, dims::Dims) where ST = fill(down(ST), dims)

Random.rand(::Type{Bit{T}}, dims::Dims) where T = Bit{T}.(rand(Bool, dims))
Random.rand(::Type{Spin{T}}, dims::Dims) where T = Spin{T}.(2 * rand(Bool, dims) .- 1)
Random.rand(::Type{Half{T}}, dims::Dims) where T = Half{T}.(rand(Bool, dims) .- 0.5)

# static array
import Base: @_inline_meta

ups(::Type{SA}) where {SA <: StaticArray} = _ups(Size(SA), SA)

@generated function _ups(::Size{s}, ::Type{SA}) where {s, SA <: StaticArray}
    T = eltype(SA)
    if T == Any
        T = Bit{Float64} # default
    end
    v = [:(up($T)) for i = 1:prod(s)]
    return quote
        @_inline_meta
        $SA(tuple($(v...)))
    end
end

downs(::Type{SA}) where {SA <: StaticArray} = _downs(Size(SA), SA)

@generated function _downs(::Size{s}, ::Type{SA}) where {s, SA <: StaticArray}
    T = eltype(SA)
    if T == Any
        T = Bit{Float64} # default
    end
    v = [:(down($T)) for i = 1:prod(s)]
    return quote
        @_inline_meta
        $SA(tuple($(v...)))
    end
end
