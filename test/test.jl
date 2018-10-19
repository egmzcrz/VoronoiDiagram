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
