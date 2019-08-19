mutable struct Beachline
    root::Union{AbstractVoronoiNode, Nothing}
    Beachline() = new(nothing)
end

mutable struct Breakpoint <: AbstractVoronoiNode
    parent::Union{AbstractVoronoiNode, Nothing}
    left::Union{AbstractVoronoiNode, Nothing}
    right::Union{AbstractVoronoiNode, Nothing}
    parabolas::NTuple{2, NTuple{2,Float64}}
    halfedge::Union{Halfedge, Nothing}
    Breakpoint(parabolas) = new(nothing, nothing, nothing,
                                         parabolas, nothing)
end

mutable struct Arc <: AbstractVoronoiNode
    parent::Union{AbstractVoronoiNode, Nothing}
    site::NTuple{2, Float64}
    circleEvent::Union{AbstractVoronoiEvent, Nothing}
    Arc(site) = new(nothing, site, nothing)
end

function getLeftBreakpoint(arc::Arc)
    currentNode = arc
    breakpoint = currentNode.parent
    while breakpoint !== nothing && breakpoint.left === currentNode
        currentNode = breakpoint
        breakpoint = currentNode.parent
    end
    return breakpoint
end

function getRightBreakpoint(arc::Arc)
    currentNode = arc
    breakpoint = currentNode.parent
    while breakpoint !== nothing && breakpoint.right === currentNode
        currentNode = breakpoint
        breakpoint = currentNode.parent
    end
    return breakpoint
end

function getLeftArc(breakpoint::Breakpoint)
    currentNode = breakpoint.left
    while isa(currentNode, Breakpoint)
        currentNode = currentNode.right
    end
    return currentNode
end

function getRightArc(breakpoint::Breakpoint)
    currentNode = breakpoint.right
    while isa(currentNode, Breakpoint)
        currentNode = currentNode.left
    end
    return currentNode
end

function setLeftChild!(breakpoint::Breakpoint, child::AbstractVoronoiNode)
    breakpoint.left = child
    child.parent = breakpoint
end

function setRightChild!(breakpoint::Breakpoint, child::AbstractVoronoiNode)
    breakpoint.right = child
    child.parent = breakpoint
end
