include("./src/VoronoiDiagram.jl");

function test()
    sites = [VoronoiDiagram.Point(rand(),rand()) for i = 1:1000]
    @time edges = VoronoiDiagram.getVoronoiDiagram(sites)

    gp = open(`gnuplot -p`, "w")
    println(gp, "
        set terminal qt noraise;
        set colorsequence podo;
        unset key; set size ratio -1;
        set tics out nomirror;
        set xrange [0:1];
        set yrange [0:1];
        plot '-' w p pt 7 lt 2, '-' w l lt 1;
        ")
    for site in sites
        println(gp, site.x, " ", site.y)
    end
    println(gp,"e")

    for e in edges
        println(gp, e.start.x, " ", e.start.y)
        println(gp, e.finish.x, " ", e.finish.y, "\n")
    end
    println(gp,"e")
    close(gp)
end
@time test()
