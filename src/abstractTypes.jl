abstract type AbstractVoronoiNode end

abstract type AbstractVoronoiEvent end

abstract type AbstractPolygonalMesh end
abstract type AbstractVertex <: AbstractPolygonalMesh end
abstract type AbstractFace <: AbstractPolygonalMesh end
abstract type AbstractHalfedge <: AbstractPolygonalMesh end
