function gnuplot(x::Int,y::Int)
    p = open(pipeline(`gnuplot`,stdout=devnull, stderr=devnull),"w")
    print(p,"set terminal epslatex standalone size $(x)cm,$(y)cm\\
            color background 'white' lw 3 font ',20'\\
            header '\\usepackage[utf8]{inputenc} \\usepackage{amsmath}';
            set output 'gptemp.tex'; set colorsequence podo;")
    return p
end
function plot(p::Base.Process)
    close(p.in)
    close(p.out)
    close(p.err)
    wait(p.closenotify)
    run(pipeline(`latex gptemp.tex`, stdout=devnull, stderr=devnull))
    run(pipeline(`dvips gptemp.dvi`, stdout=devnull, stderr=devnull))
    run(pipeline(`gs -dBATCH -dSAFER -dNOPAUSE -sDEVICE=pngalpha
            -sOutputFile="output.png" gptemp.ps`,stdout=devnull,stderr=devnull))
    run(`find -x . -name gptemp\* -delete`)
    display("image/png",read("output.png"))
end
