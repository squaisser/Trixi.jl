using OrdinaryDiffEqLowStorageRK
using OrdinaryDiffEqSDIRK
using Trixi
using Plots


###############################################################################
# semidiscretization of the relaxed Euler equations
eps = 1.0e-3
a = 1.0
equations = IncompressibleEulerRelaxationEquations1D(eps, a)
tspan = (0.0, 3.23*eps)

function initial_condition_gauss_wall(x, t, equations::IncompressibleEulerRelaxationEquations1D)
    p_eps = 0.0
    v1 = 2 * exp(-(x[1])^2 / 0.05)
    V_eps = 0.0
    return SVector(p_eps, v1, V_eps)
end

initial_condition = initial_condition_gauss_wall
source_terms = source_terms_constant

solver = DGSEM(polydeg = 4, surface_flux = flux_hll)

coordinates_min = (-1.0,)
coordinates_max = (1.0,)

mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 5,
                n_cells_max = 30_000,
                periodicity = false
                )

boundary_conditions = (;
    x_neg = boundary_condition_slip_wall,
    x_pos = boundary_condition_slip_wall
)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms,
                                    boundary_conditions = boundary_conditions
                                    )

###############################################################################
# ODE solvers, callbacks etc.

ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                            analysis_integrals=()) #NOTE: removed extra analysis errors and integrals to avoid errors

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 100,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2cons)

stepsize_callback = StepsizeCallback(cfl = 0.3)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        #save_solution,
                        stepsize_callback,  
                        )

###############################################################################
# run the simulation
steps_vis = 20
lssolver = CarpenterKennedy2N54(williamson_condition = false);
#lssolver = ImplicitEuler(autodiff=false)
sol = solve(ode, lssolver;
            dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback for explicit solvers
            ode_default_options()..., callback = callbacks, saveat = range(ode.tspan..., length=steps_vis));
println("Simulation finished with code $(sol.retcode).")

###############################################################################
# plot some results

doPlot = true
create_gif = false

if create_gif || sol.retcode == :Success
    anim = @animate for i in 1:steps_vis
        pd = PlotData1D(sol.u[i], semi)
        plot(pd["V_eps"])
    end
    mygif = gif(anim, "plot/V_eps.gif"; fps=steps_vis/3)
end

if !doPlot || sol.retcode != :Success
    println("Plotting skipped.")
else
    for i in 1:steps_vis
        pd = PlotData1D(sol.u[i], semi)
        p_v1 = plot(pd["v1"], ylim=(-2.0, 2.0))
        p_peps = plot(pd["p_eps"])
        p_Veps = plot(pd["V_eps"])
        display(plot(p_v1, p_peps, p_Veps, layout = (3, 1)))
    end
    println("Plotting finished.")
end