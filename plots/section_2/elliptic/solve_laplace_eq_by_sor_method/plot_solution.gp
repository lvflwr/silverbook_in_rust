set terminal pngcairo size 1280, 960 enhanced font ",24"

set xlabel "x"
set ylabel "y"
unset xtics
unset ytics

set pm3d map
set palette rgbformulae 21,22,23

set output "outputs/section_2/elliptic/solve_laplace_eq_by_sor_method/solution.png"
splot "outputs/section_2/elliptic/solve_laplace_eq_by_sor_method/solution.dat" u 1:2:3 notitle
