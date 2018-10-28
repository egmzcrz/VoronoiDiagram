include("maxPQ.jl")

abstract type VoronoiEvent end
abstract type VoronoiNode end

struct Point
    x::Float64
    y::Float64
end

mutable struct HalfEdge
    p::Point
    q::Point
    direction::Point
    HalfEdge(p::Point,direction::Point) = new(p,p,direction)
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
    siteL::Point
    siteR::Point
    halfedge::HalfEdge

    parent::Union{VoronoiNode,Nothing}
    childL::Union{VoronoiNode,Nothing}
    childR::Union{VoronoiNode,Nothing}

    function Breakpoint(p::Point, q::Point, halfedge::HalfEdge)
        return new(p,q,halfedge,nothing,nothing,nothing)
    end
end

mutable struct Arc <: VoronoiNode
    site::Point
    circleEvent::Union{VoronoiEvent,Nothing}
    parent::Union{VoronoiNode,Nothing}
    Arc(p::Point) = new(p,nothing,nothing)
end

# A priority queue handles events. Events have reference to points.
# They are ordered based on their point's position.
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

import Base.>
function >(e1::VoronoiEvent, e2::VoronoiEvent)
    p = e1.p
    q = e2.p
    if p.y == q.y; return p.x < q.x; end
    return p.y > q.y
end

# Parabola with focal point, f, and directrix ly evaluated at x
function parabola(x::Float64, ly::Float64, f::Point)
    a = 2*(f.y - ly)
    b = f.x*f.x + f.y*f.y - ly*ly
    return (x*x - 2*x*f.x + b)/a
end

# ly is sweepline's y-coordinate
# p and q represent the foci of the parabolas they define.
function calculateBreakpoint(p::Point, q::Point, ly::Float64)
    if p.y == q.y; return (p.x + q.x) / 2; end
    # intersection of two parabolas
    a = p.y - ly
    b = q.y - ly
    c = p.x - q.x
    d = p.y - q.y
    s1 = a*q.x - b*p.x
    s2 = sqrt((c*c + d*d)*a*b)
    
    # There are two solutions: x1, x2, we want the point whose
    # x-coordinate is in between both sites.
    x1 = (s1 - s2) / d
    x2 = (s1 + s2) / d
    return p.y < q.y ? max(x1,x2) : min(x1,x2)
end

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

function ccw(a::Point, b::Point, c::Point)
    abx = b.x - a.x
    aby = b.y - a.y
    acx = c.x - a.x
    acy = c.y - a.y

    crossProduct = abx * acy - aby * acx
    # ccw will take into account crossProduct >= 0,
    # that is, collinearity is accounted for.
    return crossProduct > -1e-15
end

function getArcOnTop(beachline::BST, site::Point)
    node = beachline.root
    while typeof(node) === Breakpoint
        x = calculateBreakpoint(node.siteL, node.siteR, site.y)
        node = site.x < x ? node.childL : node.childR
    end
    return node
end

function getLeftBreakpoint(node::Arc)
    parent = node.parent
    if parent === nothing; return nothing; end
    while parent.childL === node
        if parent.parent === nothing; return nothing; end
        node = parent
        parent = parent.parent
    end
    return parent
end
function getRightBreakpoint(node::Arc)
    parent = node.parent
    if parent === nothing; return nothing; end
    while parent.childR === node
        if parent.parent === nothing; return nothing; end
        node = parent
        parent = parent.parent
    end
    return parent
end
function getLeftArc(node::Breakpoint)
    if node === nothing; return nothing; end
    child = node.childL
    while typeof(child) === Breakpoint
        child = child.childR
    end
    return child
end
function getRightArc(node::Breakpoint)
    if node === nothing; return nothing; end
    child = node.childR
    while typeof(child) === Breakpoint
        child = child.childL
    end
    return child
end
function setRightChild!(node::Breakpoint, child::VoronoiNode)
    node.childR = child
    child.parent = node
end
function setLeftChild!(node::Breakpoint, child::VoronoiNode)
    node.childL = child
    child.parent = node
end

function checkCircleEvent(arc::Arc, events::Vector{VoronoiEvent})
    # Get arc's left and right neighbouring arcs.
    lbp = getLeftBreakpoint(arc)
    rbp = getRightBreakpoint(arc)
    if lbp === nothing || rbp === nothing; return; end
    la = getLeftArc(lbp)
    ra = getRightArc(rbp)
    if la === nothing || ra === nothing || la.site === ra.site; return; end

    # Check if edges <la,arc> and <arc,ra> will converge.
    a = la.site
    b = arc.site
    c = ra.site
    if ccw(a,b,c); return; end

    # Edges will converge at the circumcenter of the circle defined
    # by the sites of la, arc, and ra. The circle event will then take
    # place when the sweep line reaches the bottom of that circle.
    #
    # Create circle event pointing to arc.
    p = getCircleLowestPoint(a,b,c)
    circleEvent = CircleEvent(p,arc)
    # Point arc to circle event
    arc.circleEvent = circleEvent

    enqueue!(events, circleEvent)
end


function markAsFalseAlarm!(event::CircleEvent)
    event.isFalseAlarm = true
    event.arc.circleEvent = nothing
    event.arc = nothing
end

function handleSiteEvent(site::Point, beachline::BST,
                         events::Vector{VoronoiEvent}, edges::Vector{HalfEdge})
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

    # Start new edges at breakpoint (intersection between
    # current site's parabola and the arc directly on top)
    # Edges will grow in the directions parallel to the
    # bisector between site and arc's focus.
    x = site.x
    y = parabola(x, site.y, arc.site)
    start = Point(x,y)
    # A vector pointing from current site to its arc on top.
    rx = arc.site.x - site.x
    ry = arc.site.y - site.y
    # (-ry,rx) points to the left parallel to the bisector.
    # (ry,-rx) points to the right parallel to the bisector.
    el = HalfEdge(start, Point(-ry,rx))
    er = HalfEdge(start, Point(ry,-rx))
    push!(edges,el,er)


    # Replace arc-on-top with new subtree:
    #   lbp
    #  /   \
    # la   rbp
    #     /   \
    #    ma   ra
    lbp = Breakpoint(arc.site, site, el)
    rbp = Breakpoint(site, arc.site, er)
    la = Arc(arc.site)
    ma = Arc(site)
    ra = Arc(arc.site)
    setLeftChild!(lbp, la)
    setRightChild!(lbp, rbp)
    setLeftChild!(rbp, ma)
    setRightChild!(rbp, ra)
    if arc.parent !== nothing
        if arc.parent.childL === arc
            setLeftChild!(arc.parent, lbp)
        else
            setRightChild!(arc.parent, lbp)
        end
    else
        beachline.root = lbp
    end

    # Rebalancing operations?

    # Left or right arcs could eventually disappear
    # check for their circle events.
    checkCircleEvent(la,events)
    checkCircleEvent(ra,events)
end 

function handleCircleEvent(event::CircleEvent, events::Vector{VoronoiEvent},
                           edges::Vector{HalfEdge})
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

    # lbp's and rbp's edges will meet at circumcenter.
    # Finish those edges.
    x = event.p.x
    y = parabola(x, event.p.y, la.site)
    finish = Point(x,y)
    lbp.halfedge.q = finish
    rbp.halfedge.q = finish

    # lbp and rbp are the only nodes that reference the arc that will disappear.
    # One of them is the direct parent of the disappearing arc and it
    # has to be removed. The other needs to be updated to reference the new
    # growing edge between la and ra.
    innerMost = arc.parent === lbp ? rbp : lbp
    innerMost.siteL = la.site
    innerMost.siteR = ra.site
    rx = la.site.x - ra.site.x
    ry = la.site.y - ra.site.y
    edge = HalfEdge(finish, Point(-ry,rx))
    push!(edges, edge)
    innerMost.halfedge = edge

    # Save the node that is not the disappearing arc and connect it
    # to arc's grandparent.
    gparent = arc.parent.parent
    parent = arc.parent
    node = parent.childL === arc ? parent.childR : parent.childL
    if gparent.childL === parent
        setLeftChild!(gparent,node)
    else
        setRightChild!(gparent,node)
    end

    checkCircleEvent(la,events)
    checkCircleEvent(ra,events)
end

function finishEdges!(beachline::BST)
    node = beachline.root
    if typeof(node) === Arc; return; end
    breakpts::Vector{Breakpoint} = [node]
    while length(breakpts) > 0
        node = pop!(breakpts)
        childL = node.childL
        childR = node.childR
        if typeof(childL) === Breakpoint
            push!(breakpts, childL)
        end
        if typeof(childR) === Breakpoint
            push!(breakpts, childR)
        end
        e = node.halfedge
        x = e.p.x + e.direction.x 
        y = e.p.y + e.direction.y 
        e.q = Point(x,y)
    end
end

function getVoronoiDiagram(sites::Vector{Point})
    # Initialize event queue with all site events.
    events::Vector{VoronoiEvent} = []
    for site in sites
        enqueue!(events, SiteEvent(site))
    end

    # Initialize empty beachline
    beachline::BST = BST()

    # Initialize empty edge list
    edges::Vector{HalfEdge} = []

    while length(events) > 0
        event = dequeue!(events)
        if typeof(event) === SiteEvent
            handleSiteEvent(event.p, beachline, events, edges)
        elseif !event.isFalseAlarm
            handleCircleEvent(event, events, edges)
        end
    end

    finishEdges!(beachline)
    return edges
end
