include("../src/priorityQueue.jl")

# Priority queue works for any type for which
# Base.< has been overloaded.
struct Comparable
    key::Float64
    data::Int
end
import Base.<
function <(x::Comparable, y::Comparable)
    return x.key < y.key
end

function testPQ()
    arr::Vector{Comparable} = []

    for i in 1:10
        x = Comparable(rand(), 1)
        enqueue!(arr,x)
    end

    for i in eachindex(arr)
        k = dequeue!(arr)
        println(k)
    end
end

testPQ()

sites = [Point(rand(),rand()) for i = 1:10000]
@time edges = getVoronoiDiagram(sites)

p = open(`gnuplot`,"w")
println(p, "set terminal png; set output 'output.png'; unset key; set size ratio -1;
        set xrange [0:1]; set yrange [0:1]; set tics out;")
println(p, "plot '-' w p pt 6 ps .1, '-' w l")
for site in sites
    px = site.x
    py = site.y
    println(p, "$px $py")
end
println(p,"e")

for e in edges
    px = e.p.x
    py = e.p.y
    qx = e.q.x
    qy = e.q.y
    println(p, "$px $py")
    println(p, "$qx $qy\n")
end
println(p,"e")
close(p)
