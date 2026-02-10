# The same setup as tree_1d_dgsem/elixir_euler_source_terms.jl
# to verify the StructuredMesh implementation against TreeMesh

using OrdinaryDiffEqLowStorageRK
using OrdinaryDiffEq
using Trixi
#using Infiltrator
using Plots

###############################################################################
# semidiscretization of the compressible Euler equations
eps = 1e-3
a = 1.0
equations = IncompressibleEulerRelaxationEquations2D(eps, a)
tspan = (0.0, eps/sqrt(2)*0.5)

initial_condition = initial_condition_riemann
#initial_condition = initial_condition_constant
#source_terms = source_terms_homogeneous
source_terms = source_terms_constant

boundary_condition_zero_dirichlet = BoundaryConditionDirichlet((x, t, equations) -> SVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
boundary_conditions_hyperbolic = (;
                                  x_neg = boundary_condition_do_nothing,
                                  x_pos = boundary_condition_do_nothing,
                                  y_neg = boundary_condition_do_nothing,
                                  y_pos = boundary_condition_do_nothing
                                 )

solver = DGSEM(polydeg = 4, surface_flux = FluxLaxFriedrichs())

coordinates_min = (-1.0, -1.0)
coordinates_max = (1.0, 1.0)

#cells_per_dimension = (100, 100) #fails with implicit Euler due to out of memory
#cells_per_dimension = (16, 16)
#mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max,
#                periodicity = false)
initial_refinement_level = 5 #changed from 5 to 3 to reduce number of cells (out of memory)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = initial_refinement_level,
                n_cells_max = 30_000,
                periodicity = false
                )

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms,
                                    boundary_conditions = boundary_conditions_hyperbolic
                                    )

###############################################################################
# ODE solvers, callbacks etc.

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                            analysis_integrals=())

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 100,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2cons)

stepsize_callback = StepsizeCallback(cfl = 0.3)

#time_series = TimeSeriesCallback(semi, [(-0.5, -0.5), (0.5, 0.5)];
#                                 interval=5,
#                                 solution_variables=cons2cons,
#                                 filename="tseries.h5")

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback,
                        #time_series
                        )

###############################################################################
# run the simulation
steps_vis = 6
println("Starting simulation...")
#println("Initial refinement level: $(initial_refinement_level)")
lssolver = CarpenterKennedy2N54(williamson_condition = false)
#lssolver = ImplicitEuler(autodiff=false)
sol = solve(ode, lssolver;
            dt = 0.0001, # solve needs some value here but it will be overwritten by the stepsize_callback for explicit solvers
            ode_default_options()..., callback = callbacks);#, saveat = range(ode.tspan..., length=steps_vis));
println("Simulation finished with code $(sol.retcode).")

###############################################################################
# plot some results
doPlot = true
if !doPlot || sol.retcode != :Success
    println("Plotting skipped.")
else
    pd = PlotData1D(sol)
    plot!(p1, pd["p_eps"], label="2D", lw=2, legend=:best)
    plot!(p2, pd["v1"], label="2D", lw=2, legend=:best)
    plot!(p3, pd["V_eps_11"], label="2D", lw=2, legend=:best)
    display(p1)
    display(p2)
    display(p3)

    pp = plot(pd["p_eps"], title="p_eps", lw=2)
    pv1 = plot(pd["v1"], title="v1", lw=2)
    pv2 = plot(pd["v2"], title="v2", lw=2)
    pVeps_11 = plot(pd["V_eps_11"], title="V_eps_11", lw=2)
    pVeps_12 = plot(pd["V_eps_12"], title="V_eps_12", lw=2)
    pVeps_21 = plot(pd["V_eps_21"], title="V_eps_21", lw=2)
    pVeps_22 = plot(pd["V_eps_22"], title="V_eps_22", lw=2)
    display(pp)
    display(pv1)
    display(pv2)
    display(pVeps_11)
    display(pVeps_12)
    display(pVeps_21)
    display(pVeps_22)
    println("Plotting finished.")
end