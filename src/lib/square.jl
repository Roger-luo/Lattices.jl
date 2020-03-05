export Square

struct Square{BC <: Boundary, N} <: BoundedLattice{BC, 2} end

Square{BC}(n::Int) where BC = Square{BC, n}()
Square(n::Int) = Square{Periodic}(n)

Base.size(l::Square{BC, N}) where {BC, N} = (N, N)
Base.length(l::Square{BC, N}) where {BC, N} = N * N

Base.show(io::IO, x::Square{BC, N}) where {BC, N} = print(io, "Square{$BC}(", N, ")")

function to_coordinate(l::Square{BC, N}, k::Int) where {BC, N}
    @assert 1 <= k <= length(l) "site id is too large, expect 1 <= k <= $(length(l))"
    t = k - 1
    x2 = t ÷ N + 1
    t = t - (x2 - 1) * N
    x1 = t + 1
    return Coordinate(x1, x2)
end

function to_siteid(l::Square{BC, N}, x::Coordinate{2}) where {BC, N}
    x1, x2 = x
    return x1 + (x2 - 1) * N
end

Base.size(::SiteIt{<:Square{BC, N}}) where {BC, N} = (N, N)
Base.IteratorSize(::Type{<:SiteIt{Square{BC, N}}}) where {BC, N} = Base.HasShape{2}()


Base.length(::EdgeIt{<:Square{Periodic, N}}) where N = 2N * N
Base.length(::EdgeIt{<:Square{Open, N}}) where N = 2N * (N - 1)

struct SquareEdgeState{D}
    st::Int
end

Base.:(+)(st::SquareEdgeState{D}, x::Int) where D = SquareEdgeState{D}(st.st + x)

function Base.iterate(it::EdgeIt{<:Square, K}) where K
    if isodd(K)
        return iterate(it, SquareEdgeState{:-}(1))
    else
        return iterate(it, SquareEdgeState{:/}(1))
    end
end

function Base.iterate(it::EdgeIt{<:Square{Periodic, N}, K}, st::SquareEdgeState{:-}) where {K, N}
    id = st.st
    distance = K ÷ 2 + 1
    id > length(it.lattice) && return iterate(it, SquareEdgeState{:|}(1))

    x1, x2 = to_coordinate(it.lattice, id)
    y1 = x1
    y2 = mod1(x2 + distance, N)
    return (id, to_siteid(it.lattice, Coordinate(y1, y2))), st + 1
end

function Base.iterate(it::EdgeIt{<:Square{Periodic, N}, K}, st::SquareEdgeState{:|}) where {K, N}
    id = st.st
    distance = K ÷ 2 + 1
    id > length(it.lattice) && return nothing

    x1, x2 = to_coordinate(it.lattice, id)
    y1 = mod1(x1 + distance, N)
    y2 = x2
    return (id, to_siteid(it.lattice, Coordinate(y1, y2))), st + 1
end

function Base.iterate(it::EdgeIt{<:Square{Periodic, N}, K}, st::SquareEdgeState{:/}) where {K, N}
    id = st.st
    distance = K ÷ 2
    id > length(it.lattice) && return iterate(it, SquareEdgeState{:\}(1))

    x1, x2 = to_coordinate(it.lattice, id)
    y1 = mod1(x1 + distance, N)
    y2 = mod1(x2 + distance, N)
    return (id, to_siteid(it.lattice, Coordinate(y1, y2))), st + 1
end

function Base.iterate(it::EdgeIt{<:Square{Periodic, N}, K}, st::SquareEdgeState{:\}) where {K, N}
    id = st.st
    distance = K ÷ 2
    id > length(it.lattice) && return nothing

    x1, x2 = to_coordinate(it.lattice, id)
    y1 = mod1(x1 - distance, N)
    y2 = mod1(x2 + distance, N)
    return (id, to_siteid(it.lattice, Coordinate(y1, y2))), st + 1
end

## Open Boundary
function Base.iterate(it::EdgeIt{<:Square{Open, N}, K}, st::SquareEdgeState{:-}) where {K, N}
    id = st.st
    distance = K ÷ 2 + 1

    x1, x2 = to_coordinate(it.lattice, id)
    y1 = x1
    y2 = x2 + distance

    y2 > N && return iterate(it, SquareEdgeState{:|}(1))
    return (id, to_siteid(it.lattice, Coordinate(y1, y2))), st + 1
end

function Base.iterate(it::EdgeIt{<:Square{Open, N}, K}, st::SquareEdgeState{:|}) where {K, N}
    id = st.st
    distance = K ÷ 2 + 1
    id > length(it.lattice) && return nothing

    x1, x2 = to_coordinate(it.lattice, id)
    y1 = x1 + distance
    y2 = x2

    y1 > N && return (id, to_siteid(it.lattice, Coordinate(y1, y2))), st + distance + 1
    return (id, to_siteid(it.lattice, Coordinate(y1, y2))), st + 1
end

function Base.iterate(it::EdgeIt{<:Square{Open, N}, K}, st::SquareEdgeState{:/}) where {K, N}
    id = st.st
    distance = K ÷ 2

    x1, x2 = to_coordinate(it.lattice, id)
    y1 = x1 + distance
    y2 = x2 + distance

    y1 > N && return (id, to_siteid(it.lattice, Coordinate(y1, y2))), st + distance + 1
    y2 > N && return iterate(it, SquareEdgeState{:\}(2))
    return (id, to_siteid(it.lattice, Coordinate(y1, y2))), st + 1
end

function Base.iterate(it::EdgeIt{<:Square{Open, N}, K}, st::SquareEdgeState{:\}) where {K, N}
    id = st.st
    distance = K ÷ 2
 
    x1, x2 = to_coordinate(it.lattice, id)
    y1 = x1 - distance
    y2 = x2 + distance

    x1 == N && return (id, to_siteid(it.lattice, Coordinate(y1, y2))), st + distance + 1 # skip
    x2 > N && return nothing
    return (id, to_siteid(it.lattice, Coordinate(y1, y2))), st + 1
end

# Face Iterator
Base.length(::FaceIt{<:Square{Periodic, N}, K}) where {N, K} = N * N
Base.size(::FaceIt{<:Square{BC, N}}) where {BC, N} = (N, N)
Base.IteratorSize(::Type{<:FaceIt{Square{BC, N}}}) where {BC, N} = Base.HasShape{2}()


function Base.iterate(it::FaceIt{<:Square, K}) where K
    if isodd(K)
        return iterate(it, SquareEdgeState{:-}(1))
    else
        return iterate(it, SquareEdgeState{:/}(1))
    end
end

function Base.iterate(it::FaceIt{<:Square{Periodic, N}, K}, st::SquareEdgeState{:-}) where {N, K}
    id = st.st
    id > N * N && return nothing
    distance = K ÷ 2 + 1
    a1, a2 = to_coordinate(it.lattice, id)

    b1, b2 = mod1(a1 + distance, N), a2
    c1, c2 = mod1(a1 + distance, N), mod1(a2 + distance, N)
    d1, d2 = a1, mod1(a2 + distance, N)

    (id,
        to_siteid(it.lattice, Coordinate(b1, b2)),
        to_siteid(it.lattice, Coordinate(c1, c2)),
        to_siteid(it.lattice, Coordinate(d1, d2)),
    ), st + 1
end


function Base.iterate(it::FaceIt{<:Square{Periodic, N}, K}, st::SquareEdgeState{:/}) where {N, K}
    id = st.st
    id > N * N && return nothing
    distance = K ÷ 2
    a1, a2 = to_coordinate(it.lattice, id)

    b1, b2 = mod1(a1 + distance, N), mod1(a2 + distance, N)
    c1, c2 = mod1(a1 + 2distance, N), a2
    d1, d2 = mod1(a1 + distance, N), mod1(a2 - distance, N)

    (id,
        to_siteid(it.lattice, Coordinate(b1, b2)),
        to_siteid(it.lattice, Coordinate(c1, c2)),
        to_siteid(it.lattice, Coordinate(d1, d2)),
    ), st + 1
end
