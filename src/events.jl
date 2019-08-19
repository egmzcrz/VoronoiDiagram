import Base: >, <

struct SiteEvent <: AbstractVoronoiEvent
    p::NTuple{2, Float64}
end

mutable struct CircleEvent <: AbstractVoronoiEvent
    p::NTuple{2, Float64}
    arc::Union{AbstractVoronoiNode, Nothing}
    CircleEvent(p, arc) = new(p, arc)
end

function >(e1::AbstractVoronoiEvent, e2::AbstractVoronoiEvent)
    p = e1.p
    q = e2.p
    if p[2] == q[2]; return p[1] < q[1]; end
    return p[2] > q[2]
end

# TODO: check this function.
function <(e1::AbstractVoronoiEvent, e2::AbstractVoronoiEvent)
    p = e1.p
    q = e2.p
    if p[2] == q[2]; return p[1] < q[1]; end
    return p[2] > q[2]
end

function enqueue!(heap::Vector{T}, key::T) where {T}
    push!(heap,key)
    swim!(heap)
end

function dequeue!(heap::Vector{T}) where {T}
    if length(heap) > 1
        last = length(heap)
        heap[1], heap[last] = heap[last], heap[1]
    end
    removedKey = pop!(heap)
    sink!(heap)
    return removedKey
end

up(i::Int) = i >>> 1
left(i::Int) = i << 1
right(i::Int) = (i << 1) + 1

function swim!(heap::Vector{T}) where {T}
    node = length(heap)
    upper = up(node)
    while node > 1 && heap[node] < heap[upper]
        heap[node], heap[upper] = heap[upper], heap[node]
        node = upper
        upper = up(node)
    end
end

function sink!(heap::Vector{T}) where {T}
    node = 1
    lft = left(node)
    rght = right(node)
    last = length(heap)
    while lft <= last
        minChild = lft
        if lft < last && heap[rght] < heap[lft]; minChild = rght end
        if heap[node] < heap[minChild]; break; end
        heap[node], heap[minChild] = heap[minChild], heap[node]
        node = minChild
        lft = left(node)
        rght = right(node)
    end
end
