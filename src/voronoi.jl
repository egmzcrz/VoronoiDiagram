include("abstractTypes.jl")
include("diagram.jl")
include("beachline.jl")
include("events.jl")

const Point = NTuple{2, Float64}

# Parabola with focal point f, and directrix ly evaluated at x
function parabola(x::Float64, ly::Float64, f::Point)
    fx = f[1]; fy = f[2]
    a = fx - x
    b = fy + ly
    c = fy - ly
    return (b + a * a / c) / 2
end

# fa and fb are focal points, ly is the directrix.
function intersectParabolas(fa::Point, fb::Point, ly::Float64)
    fax = fa[1]; fay = fa[2]
    fbx = fb[1]; fby = fb[2]

    if fay == fby; return fax < fbx ? (fax + fbx) / 2 : Inf; end
    # intersection of two parabolas
    a = fay - ly
    b = fby - ly
    c = fax - fbx
    d = fay - fby
    s1 = a*fbx - b*fax
    s2 = sqrt((c*c + d*d)*a*b)

    # There are two solutions: x1 = (s1-s2)/d and x2 = (s1+s2)/d.
    # Also, x1 < x2 iff d > 0.
    # So x1 will throw the point to the left or to the right of the
    # lowest parabola, which is exactly what we want.
    return (s1 - s2) / d
end

# a, b and c define a circumcircle. It's center is a Voronoi vertex.
function getCircleBottomPoint(a::Point, b::Point, c::Point)
    ax = a[1]; ay = a[2]
    bx = b[1]; by = b[2]
    cx = c[1]; cy = c[2]
    ux = bx - ax; uy = by - ay
    vx = cx - bx; vy = cy - by
    wx = cx - ax; wy = cy - ay

    vw = vx*wx + vy*wy
    uxw = ux*wy - uy*wx
    s = vw/uxw

    # Circumcenter
    x = (ax + bx - s*uy)/2
    y = (ay + by + s*ux)/2
    # Radius
    rx = ax - x
    ry = ay - y
    r = sqrt(rx*rx + ry*ry)
    return (x, y-r)
end

# ccw checks if the orientation of the ordered triplet (a,b,c)
# is counter-clockwise.
# Collinearity is when v×w = 0, so numerical errors should be considered.
function ccw(a::Point, b::Point, c::Point)
    ax = a[1]; ay = a[2]
    bx = b[1]; by = b[2]
    cx = c[1]; cy = c[2]
    abx = bx - ax
    aby = by - ay
    acx = cx - ax
    acy = cy - ay

    crossProduct = abx * acy - aby * acx
    return crossProduct > -1e-15
end

# Search for the arc directly on top of current site.
function getSplittingArc(beachline::Beachline, site::Point)
    currentNode = beachline.root
    while isa(currentNode, Breakpoint)
        parabolas = currentNode.parabolas
        x = intersectParabolas(parabolas[1], parabolas[2], site[2])
        currentNode = site[1] < x ? currentNode.left : currentNode.right
    end
    return currentNode
end

function checkCircleEvent(la::Arc, ma::Arc, ra::Arc,
                          events::Vector{AbstractVoronoiEvent})
    # Check if edges <la,ma> and <ma,ra> converge.
    a = la.site
    b = ma.site
    c = ra.site
    if ccw(a,b,c); return; end

    # Edges will converge at the circumcenter defined by a,b,c.
    # The circle event will then take place when the sweep line
    # reaches the bottom of that circle.
    p = getCircleBottomPoint(a, b, c)
    # Create circle event pointing to arc.
    circleEvent = CircleEvent(p, ma)
    # Point arc to circle event
    ma.circleEvent = circleEvent

    enqueue!(events, circleEvent)
end

function markAsFalseAlarm!(event::CircleEvent)
    event.arc.circleEvent = nothing
    event.arc = nothing
end

function handleSiteEvent!(beachline::Beachline,
                          events::Vector{AbstractVoronoiEvent},
                          diagram::PolygonalMesh,
                          site::Point)
    face = Face(site)
    push!(diagram.faces, face)

    if beachline.root === nothing
        beachline.root = Arc(site)
        return
    end

    # Get arc vertically above site (this arc will split in two) and check if
    # it points to a circle event. If it does then this event is a false alarm.
    arc = getSplittingArc(beachline, site)
    if arc.circleEvent !== nothing
        markAsFalseAlarm!(arc.circleEvent)
    end

    # Replace arc with following subtree:
    #   lbp ~~~~> halfedge
    #  /   \
    # la   rbp ~~~~> twin-halfedge
    #     /   \
    #    ma   ra
    lbp = Breakpoint((arc.site, site))
    rbp = Breakpoint((site, arc.site))
    la = Arc(arc.site)
    ma = Arc(site)
    ra = Arc(arc.site)
    setLeftChild!(lbp, la)
    setRightChild!(lbp, rbp)
    setLeftChild!(rbp, ma)
    setRightChild!(rbp, ra)
    if arc !== beachline.root
        v = arc.parent
        if v.left === arc
            setLeftChild!(v, lbp)
        else
            setRightChild!(v, lbp)
        end
    else
        beachline.root = lbp
    end

    # TODO: Rebalancing operations.

    # Create halfedge records in diagram to be traced by the new breakpoints.
    halfedge = Halfedge()
    twin = Halfedge()
    halfedge.twin = twin
    twin.twin = halfedge
    lbp.halfedge = halfedge
    rbp.halfedge = twin

    # Create a new face and link it accordingly.
    face.outerComponent = halfedge
    halfedge.incidentFace = face

    # la or ra could eventually disappear, check for potential circle events.
    llbp = getLeftBreakpoint(la)
    rrbp = getRightBreakpoint(ra)
    if llbp !== nothing
        lla = getLeftArc(llbp)
        if lla.site !== ma.site; checkCircleEvent(lla, la, ma, events); end
    end
    if rrbp !== nothing
        rra = getRightArc(rrbp)
        if ma.site !== rra.site; checkCircleEvent(ma, ra, rra, events); end
    end
end

function handleCircleEvent!(events::Vector{AbstractVoronoiEvent},
                           diagram::PolygonalMesh,
                           arc::Arc)

    # Delete disappearing arc and update breakpoints.
    #
    # Find its left and right neighbouring arcs.
    lbp = getLeftBreakpoint(arc)
    rbp = getRightBreakpoint(arc)
    la = getLeftArc(lbp)
    ra = getRightArc(rbp)
    #
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
    dp = arc.parent # direct parent
    gp = dp.parent # grandparent
    #
    # Update inner-most node's breakpoint to <la, ra>.
    innerMost = dp === rbp ? lbp : rbp
    innerMost.parabolas = (la.site, ra.site)
    #
    # Update grandparent node to point to la or ra, accordingly.
    # (tree example above shows it should point to ra, but in general we
    # should check which is the case to handle).
    node = dp.left === arc ? dp.right : dp.left
    gp.left === dp ? setLeftChild!(gp, node) : setRightChild!(gp, node)
    #
    # TODO: Rebalancing operations.
    #
    # Mark la and ra as false alarms.
    if la.circleEvent !== nothing
        markAsFalseAlarm!(la.circleEvent)
    end
    if ra.circleEvent !== nothing
        markAsFalseAlarm!(ra.circleEvent)
    end

    # Add the center of the circle as a Voronoi vertex.
    #
    # Get the Voronoi vertex.
    p = arc.circleEvent.p
    x = p[1]; ly = p[2]
    y = parabola(x, ly, arc.site) # la's, ra's or arc's site will do.
    vv = Vertex((x,y))
    #
    # Get halfedges involved.
    lhe = lbp.halfedge
    rhe = rbp.halfedge
    mhe = Halfedge() # halfedge <la, ra>.
    twin = Halfedge() # mhe's twin.
    mhe.twin = twin
    twin.twin = mhe
    #
    # Link vertex and halfedges accordingly.
    vv.incidentHalfedge = mhe

    lhe.twin.origin = vv
    rhe.twin.origin = vv
    mhe.origin = vv

    lhe.next = rhe.twin
    rhe.twin.prev = lhe

    rhe.next = mhe
    mhe.prev = rhe

    twin.next = lhe.twin
    lhe.twin.prev = twin

    innerMost.halfedge = mhe

    # Push finished halfedges
    push!(diagram.halfedges, lhe, rhe, mhe.twin)
    push!(diagram.vertices, vv)

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

function finishEdges!(beachline::Beachline, diagram::PolygonalMesh, boundingBox::Vector{Halfedge})
    node = beachline.root
    if isa(node, Arc); return; end
    breakpts::Vector{Breakpoint} = [node]
    #
    # Traverse remaining nodes of beachline (infinite halfedges).
    while length(breakpts) > 0
        breakpoint = pop!(breakpts)

        halfedge = breakpoint.halfedge
        #
        # Halfedge's direction (clockwise rotation), s.
        parabolas = breakpoint.parabolas
        lsite = parabolas[1]
        rsite = parabolas[2]
        nx = rsite[1] - lsite[1]
        ny = rsite[2] - lsite[2]
        s = (ny, -nx)

        v = Vertex(halfedge.origin.coordinates .+ 100 .* s)
        v.incidentHalfedge = halfedge.twin
        halfedge.twin.origin = v
        tmp = halfedge
        while tmp.prev !== nothing
            tmp = tmp.prev
        end
        halfedge.next = tmp
        tmp.prev = halfedge

        # Push finished halfedges
        push!(diagram.halfedges, halfedge)

        lchild = breakpoint.left
        rchild = breakpoint.right
        if isa(lchild, Breakpoint)
            push!(breakpts, lchild)
        end
        if isa(rchild, Breakpoint)
            push!(breakpts, rchild)
        end
    end
end


function boundingBox(pts::Vector{Point})
    n::Int = length(pts)
    bbox::Vector{Halfedge} = [Halfedge() for i in 1:n]
    for i in 1:n
        halfedge = bbox[i]
        halfedge.origin = Vertex(pts[i])
        halfedge.next = bbox[mod(i, n)+ 1]
        halfedge.prev = bbox[mod(n-1 + i-1, n) + 1]
        # TODO: add twin information
    end
    return bbox
end


function getVoronoiDiagram(sites::Vector{Point})
    # Initialize event queue with all site events.
    events::Vector{AbstractVoronoiEvent} = []
    for site in sites
        enqueue!(events, SiteEvent(site))
    end

    # Initialize empty beachline.
    beachline::Beachline = Beachline()

    # Initialize empty Voronoi diagram.
    diagram::PolygonalMesh = PolygonalMesh()

    while length(events) > 0
        e::AbstractVoronoiEvent = dequeue!(events)
        if isa(e, SiteEvent)
            handleSiteEvent!(beachline, events, diagram, e.p)
        else # isa circle event
            # Check if circle event is a false alarm (doesn't point to an arc).
            arc = e.arc
            if arc !== nothing
                handleCircleEvent!(events, diagram, arc)
            end
        end
    end


    diagram.faces[1].outerComponent = diagram.faces[2].outerComponent.twin
    diagram.faces[1].outerComponent.incidentFace = diagram.faces[1]


    # Make edges fit inside bounding box.
    bbox = boundingBox([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
    finishEdges!(beachline, diagram, bbox)

    # Connect every halfedge to a face.
    for face in diagram.faces
        halfedge = face.outerComponent
        while true
            halfedge.incidentFace = face

            halfedge = halfedge.next
            if halfedge === face.outerComponent
                break
            end
        end
    end

    return diagram
end
