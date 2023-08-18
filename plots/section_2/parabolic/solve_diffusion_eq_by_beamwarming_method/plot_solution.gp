set terminal pngcairo size 1280, 960 enhanced font ",24"

set xlabel "x"
set ylabel "u"

set output "outputs/section_2/parabolic/solve_diffusion_eq_by_beamwarming_method/solution.png"
plot [-1:1] for [i=0:*] "outputs/section_2/parabolic/solve_diffusion_eq_by_beamwarming_method/solution.dat" index i u 2:3 w l lw 3 title columnhead(1)
