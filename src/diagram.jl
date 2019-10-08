mutable struct Vertex <: AbstractVertex
    coordinates::NTuple{2, Float64} # Coordinates (x,y).
    incidentHalfedge::Union{AbstractHalfedge, Nothing} # Halfedge starting at this vertex.
    Vertex(p) = new(p, nothing)
end

mutable struct Halfedge <: AbstractHalfedge
    origin::Union{AbstractVertex, Nothing} # Its origin.
    twin::Union{AbstractHalfedge, Nothing} # Opposite halfedge.
    next::Union{AbstractHalfedge, Nothing} # Next halfedge ccw
    prev::Union{AbstractHalfedge, Nothing} # Next halfedge cw.
    incidentFace::Union{AbstractFace, Nothing} # Face it bounds.
    Halfedge() = new(nothing, nothing, nothing, nothing, nothing)
end

mutable struct Face <: AbstractFace
    site::NTuple{2, Float64}
    outerComponent::Union{AbstractHalfedge, Nothing}
    area::Float64
    Face(site) = new(site, nothing, 0.0)
end

struct PolygonalMesh <: AbstractPolygonalMesh
    vertices::Vector{Vertex}
    halfedges::Vector{Halfedge}
    faces::Vector{Face}
    PolygonalMesh() = new([], [], [])
end
