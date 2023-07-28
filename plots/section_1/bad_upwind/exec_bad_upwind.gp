set terminal pngcairo size 1280, 960 enhanced font ",24"

set xlabel "x"
set ylabel "u"

set output "outputs/section_1/bad_upwind/exec_bad_upwind.png"
plot [-1:1] for [i=0:*] "outputs/section_1/bad_upwind/exec_bad_upwind.dat" index i u 2:3 w l lw 3 title columnhead(1)
