module VoronoiDiagram

export getVoronoiDiagram

include("maxPQ.jl")

abstract type VoronoiEvent end
abstract type VoronoiNode end

struct Point
    x::Float64
    y::Float64
end

# A half-edge grows parallel to the bisector
# between two sites. If half-edge grows in the direction
# v, then its opposite grows in th direction -v
mutable struct Halfedge
    start::Point
    finish::Point
    opposite::Halfedge
    Halfedge(p) = new(p,p)
end

# The beachline is represented by a binary search tree, it has
# inner-nodes (breakpoints) and leaf-nodes (arcs).
# - Inner nodes keep track of breakpoints. They also point to growing
# half-edges.
# - Leaf nodes keep track of arcs of parabola. They also point to
# future circle events.
mutable struct BST
    root::Union{VoronoiNode,Nothing}
    BST() = new(nothing)
end

mutable struct Breakpoint <: VoronoiNode
    lsite::Point
    rsite::Point
    halfedge::Halfedge
    parent::Union{VoronoiNode,Nothing}
    lchild::Union{VoronoiNode,Nothing}
    rchild::Union{VoronoiNode,Nothing}
    Breakpoint(lsite, rsite, halfedge) = new(lsite, rsite, halfedge,
                                             nothing,nothing,nothing)
end

mutable struct Arc <: VoronoiNode
    site::Point
    circleEvent::Union{VoronoiEvent,Nothing}
    parent::Union{VoronoiNode,Nothing}
    Arc(site) = new(site,nothing,nothing)
end

# A priority queue handles the events. Events have references to points.
# They are ordered based on their point's location.
# At site events arcs are created.
# At circle events arcs vanish.
struct SiteEvent <: VoronoiEvent
    p::Point
end

mutable struct CircleEvent <: VoronoiEvent
    p::Point
    arc::Union{VoronoiNode,Nothing}
    isFalseAlarm::Bool
    CircleEvent(p::Point, arc::VoronoiNode) = new(p,arc,false)
end

# Overload ">" to allow comparing Voronoi events.
# This is needed for the priority queue to work with our types.
import Base.>
function >(e1::VoronoiEvent, e2::VoronoiEvent)
    p = e1.p
    q = e2.p
    if p.y == q.y; return p.x < q.x; end
    return p.y > q.y
end

# Parabola with focal point f, and directrix ly evaluated at x
function parabola(x::Float64, ly::Float64, f::Point)
    a::Float64 = f.x - x
    return 0.5 * (f.y + ly + a * a / (f.y - ly))
end

# fa and fb are focal points, ly is the directrix.
function getBreakpoint(fa::Point, fb::Point, ly::Float64)
    if fa.y == fb.y; return (fa.x + fb.x) / 2; end
    # intersection of two parabolas
    a = fa.y - ly
    b = fb.y - ly
    c = fa.x - fb.x
    d = fa.y - fb.y
    s1 = a*fb.x - b*fa.x
    s2 = sqrt((c*c + d*d)*a*b)
    
    # There are two solutions: x1 = (s1-s2)/d and x2 = (s1+s2)/d.
    # Also, x1 < x2 iff d > 0.
    # So x1 will throw the point to the left or to the right of the
    # lowest parabola, which is exactly what we want.
    return (s1 - s2) / d
end

# a, b and c define a circumcircle.
function getCircleLowestPoint(a::Point, b::Point, c::Point)
    ux = b.x - a.x; uy = b.y - a.y
    vx = c.x - b.x; vy = c.y - b.y
    wx = c.x - a.x; wy = c.y - a.y

    vw = vx*wx + vy*wy
    uxw = ux*wy - uy*wx
    s = vw/uxw

    # Circumcenter
    x = (a.x + b.x - s*uy)/2
    y = (a.y + b.y + s*ux)/2
    # Radius
    rx = a.x - x
    ry = a.y - y
    r = sqrt(rx*rx + ry*ry)
    return Point(x,y-r)
end

# ccw checks if the orientation of the ordered triplet (a,b,c)
# is counter-clockwise.
# Collinearity is when v×w = 0, so numerical errors should be considered.
function ccw(a::Point, b::Point, c::Point)
    abx = b.x - a.x
    aby = b.y - a.y
    acx = c.x - a.x
    acy = c.y - a.y

    crossProduct = abx * acy - aby * acx
    return crossProduct > -1e-15
end

# Search for the arc directly on top of current site.
function getArcOnTop(beachline::BST, site::Point)
    node = beachline.root
    while typeof(node) === Breakpoint
        x = getBreakpoint(node.lsite, node.rsite, site.y)
        node = site.x < x ? node.lchild : node.rchild
    end
    return node
end

function checkCircleEvent(la::Arc, ma::Arc, ra::Arc,
                          events::Vector{VoronoiEvent})
    # Check if edges <la,ma> and <ma,ra> converge.
    a = la.site
    b = ma.site
    c = ra.site
    if ccw(a,b,c); return; end

    # Edges will converge at the circumcenter defined by a,b,c.
    # The circle event will then take place when the sweep line
    # reaches the bottom of that circle.
    p = getCircleLowestPoint(a,b,c)
    # Create circle event pointing to arc.
    circleEvent = CircleEvent(p,ma)
    # Point arc to circle event
    ma.circleEvent = circleEvent

    enqueue!(events, circleEvent)
end

function markAsFalseAlarm!(event::CircleEvent)
    event.isFalseAlarm = true
    event.arc.circleEvent = nothing
    event.arc = nothing
end

function handleSiteEvent!(beachline::BST,
                          events::Vector{VoronoiEvent},
                          event::SiteEvent)
    site = event.p
    if beachline.root === nothing
        beachline.root = Arc(site)
        return
    end

    # Get arc directly on top of site and check if
    # it points to a circle event. If it does
    # then this event is a false alarm.
    arc = getArcOnTop(beachline, site)
    if arc.circleEvent !== nothing
        markAsFalseAlarm!(arc.circleEvent)
    end

    # Start new half-edges at breakpoint (intersection between
    # current site's parabola and the arc directly on top)
    # Edges will grow in the directions parallel to the
    # bisector between site and arc's site.
    x = site.x
    y = parabola(x, site.y, arc.site)
    start = Point(x,y)
    lhe = Halfedge(start)
    rhe = Halfedge(start)
    lhe.opposite = rhe
    rhe.opposite = lhe


    # Replace arc with following subtree:
    #   lbp ~~~~> lhe
    #  /   \
    # la   rbp ~~~~> rhe
    #     /   \
    #    ma   ra
    lbp = Breakpoint(arc.site, site, lhe)
    rbp = Breakpoint(site, arc.site, rhe)
    la = Arc(arc.site)
    ma = Arc(site)
    ra = Arc(arc.site)
    setLeftChild!(lbp, la)
    setRightChild!(lbp, rbp)
    setLeftChild!(rbp, ma)
    setRightChild!(rbp, ra)
    if arc.parent !== nothing
        if arc.parent.lchild === arc
            setLeftChild!(arc.parent, lbp)
        else
            setRightChild!(arc.parent, lbp)
        end
    else
        beachline.root = lbp
    end

    # TODO: Rebalancing operations.

    # Left or right arcs could eventually disappear,
    # check for potential circle events.
    llbp = getLeftBreakpoint(la)
    rrbp = getRightBreakpoint(ra)
    if llbp !== nothing
        lla = getLeftArc(llbp)
        if lla.site !== ma.site; checkCircleEvent(lla,la,ma,events); end
    end
    if rrbp !== nothing
        rra = getRightArc(rrbp)
        if ma.site !== rra.site; checkCircleEvent(ma,ra,rra,events); end
    end
end 

# Finishes edges that already start at a vertex or starts them at one.
# An edge that starts at a voronoi vertex shouldn't have an opposite
# edge since a vertex is like a source of edges.
function connectEdge2Vertex!(edge::Halfedge, p::Point, edges::Vector{Halfedge})
    if edge.opposite !== edge
        edge.start = p
        edge.opposite = edge
    else
        edge.finish = p
        push!(edges, edge)
    end
end

function handleCircleEvent!(edges::Vector{Halfedge},
                           events::Vector{VoronoiEvent},
                           event::CircleEvent)
    if event.isFalseAlarm; return; end
    # The arc that will disappear.
    arc = event.arc

    # Find its left and right neighbouring arcs
    lbp = getLeftBreakpoint(arc)
    rbp = getRightBreakpoint(arc)
    la = getLeftArc(lbp)
    ra = getRightArc(rbp)

    if la.circleEvent !== nothing
        markAsFalseAlarm!(la.circleEvent)
    end
    if ra.circleEvent !== nothing
        markAsFalseAlarm!(ra.circleEvent)
    end

    # Get the point of intersection.
    x = event.p.x
    ly = event.p.y
    y = parabola(x, ly, la.site) # la's, ra's or arc's site will do.
    circumcenter = Point(x,y)

    # lbp and rbp's breakpoints will be updated, so their edges will stop
    # existing. But the opposite edges will keep on growing.
    le = lbp.halfedge.opposite
    re = rbp.halfedge.opposite
    connectEdge2Vertex!(le, circumcenter, edges)
    connectEdge2Vertex!(re, circumcenter, edges)

    # lbp and rbp will be updated. One of them is the direct parent of arc.
    # The one that is the direct parent can be removed. The other one
    # needs to be updated to manage the breakpoint <la,ra>.
    #     lbp: <la,arc>            lbp: <la,ra>
    #    /    \                   /   \
    #   x     gp                 x    gp
    #  / \    / \      ==>      / \   / \
    # x  la  rbp x             x  la ra  x
    #       /  \
    #     arc   ra
    innerMost = arc.parent === rbp ? lbp : rbp
    innerMost.lsite = la.site
    innerMost.rsite = ra.site
    gp = arc.parent.parent # grandparent
    dp = arc.parent # direct parent
    node = dp.lchild === arc ? dp.rchild : dp.lchild
    gp.lchild === dp ? setLeftChild!(gp,node) : setRightChild!(gp,node)
    newEdge = Halfedge(circumcenter)
    newEdge.opposite = newEdge # newEdge starts at a voronoi vertex.
    innerMost.halfedge = newEdge

    # Check for potential circle events.
    llbp = getLeftBreakpoint(la)
    rrbp = getRightBreakpoint(ra)
    if llbp !== nothing
        lla = getLeftArc(llbp)
        checkCircleEvent(lla,la,ra,events);
    end
    if rrbp !== nothing
        rra = getRightArc(rrbp)
        checkCircleEvent(la,ra,rra,events);
    end
end

function clipEdges!(edges::Vector{Halfedge}, bbox::Vector{NTuple{2,Point}})
end

function finishEdges!(beachline::BST, edges::Vector{Halfedge})
    node = beachline.root
    if typeof(node) === Arc; return; end
    breakpts::Vector{Breakpoint} = [node]
    while length(breakpts) > 0
        node = pop!(breakpts)
        lsite = node.lsite
        rsite = node.rsite
        # Direction vector.
        rx = rsite.x - lsite.x
        ry = rsite.y - lsite.y
        # Grow half-edge parallel to direction vector.
        edge = node.halfedge
        x = edge.start.x + 5*ry
        y = edge.start.y - 5*rx
        edge.finish = Point(x,y)
        push!(edges, edge)

        lchild = node.lchild
        rchild = node.rchild
        if typeof(lchild) === Breakpoint
            push!(breakpts, lchild)
        end
        if typeof(rchild) === Breakpoint
            push!(breakpts, rchild)
        end
    end
end

function getVoronoiDiagram(sites::Vector{Point})
    # Initialize event queue with all site events.
    events::Vector{VoronoiEvent} = []
    for site in sites
        enqueue!(events, SiteEvent(site))
    end

    # Initialize empty beachline.
    beachline::BST = BST()

    # Initialize empty edge list
    edges::Vector{Halfedge} = []

    while length(events) > 0
        event::VoronoiEvent = dequeue!(events)
        if typeof(event) === SiteEvent
            handleSiteEvent!(beachline, events, event)
        else
            handleCircleEvent!(edges, events, event)
        end
    end

    # Make edges fit inside bounding box.
    #clipEdges!(edges,bbox)
    finishEdges!(beachline, edges)
    return edges
end

###################################################################
# BST utility functions
###################################################################
function getLeftBreakpoint(node::Arc)
    parent = node.parent
    if parent === nothing; return nothing; end
    while parent.lchild === node
        if parent.parent === nothing; return nothing; end
        node = parent
        parent = parent.parent
    end
    return parent
end
function getRightBreakpoint(node::Arc)
    parent = node.parent
    if parent === nothing; return nothing; end
    while parent.rchild === node
        if parent.parent === nothing; return nothing; end
        node = parent
        parent = parent.parent
    end
    return parent
end
function getLeftArc(node::Breakpoint)
    if node === nothing; return nothing; end
    child = node.lchild
    while typeof(child) === Breakpoint
        child = child.rchild
    end
    return child
end
function getRightArc(node::Breakpoint)
    if node === nothing; return nothing; end
    child = node.rchild
    while typeof(child) === Breakpoint
        child = child.lchild
    end
    return child
end
function setRightChild!(node::Breakpoint, child::VoronoiNode)
    node.rchild = child
    child.parent = node
end
function setLeftChild!(node::Breakpoint, child::VoronoiNode)
    node.lchild = child
    child.parent = node
end
###################################################################

end # module
