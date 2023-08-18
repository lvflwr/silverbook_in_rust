var sourcesIndex = JSON.parse('{\
"bad_upwind":["",[],["input.rs","lib.rs","output.rs","upwind_solver.rs"]],\
"elliptic":["",[["solver",[],["point_jacobi_solver.rs","sor_solver.rs"]]],["input.rs","lib.rs","output.rs","solver.rs"]],\
"linear_hyperbolic":["",[["math",[],["trinomial_eq.rs"]],["solver",[],["beamwarming_solver.rs","ftcs_solver.rs","lax_solver.rs","laxwendroff_solver.rs","leapfrog_solver.rs","maccormack_solver.rs","upwind_solver.rs"]]],["input.rs","lib.rs","math.rs","output.rs","solver.rs"]],\
"parabolic":["",[["math",[],["trinomial_eq.rs"]],["solver",[],["beamwarming_solver.rs","ftcs_solver.rs"]]],["input.rs","lib.rs","math.rs","output.rs","solver.rs"]],\
"solve_diffusion_eq_by_beamwarming_method":["",[],["solve_diffusion_eq_by_beamwarming_method.rs"]],\
"solve_diffusion_eq_by_ftcs_method":["",[],["solve_diffusion_eq_by_ftcs_method.rs"]],\
"solve_laplace_eq_by_point_jacobi_method":["",[],["solve_laplace_eq_by_point_jacobi_method.rs"]],\
"solve_laplace_eq_by_sor_method":["",[],["solve_laplace_eq_by_sor_method.rs"]],\
"solve_transport_eq_by_bad_upwind_method":["",[],["solve_transport_eq_by_bad_upwind_method.rs"]],\
"solve_transport_eq_by_good_upwind_method":["",[],["solve_transport_eq_by_good_upwind_method.rs"]],\
"solve_wave_eq_by_beamwarming_method":["",[],["solve_wave_eq_by_beamwarming_method.rs"]],\
"solve_wave_eq_by_ftcs_method":["",[],["solve_wave_eq_by_ftcs_method.rs"]],\
"solve_wave_eq_by_lax_method":["",[],["solve_wave_eq_by_lax_method.rs"]],\
"solve_wave_eq_by_laxwendroff_method":["",[],["solve_wave_eq_by_laxwendroff_method.rs"]],\
"solve_wave_eq_by_leapfrog_method":["",[],["solve_wave_eq_by_leapfrog_method.rs"]],\
"solve_wave_eq_by_maccormack_method":["",[],["solve_wave_eq_by_maccormack_method.rs"]],\
"solve_wave_eq_by_upwind_method":["",[],["solve_wave_eq_by_upwind_method.rs"]]\
}');
createSourceSidebar();
