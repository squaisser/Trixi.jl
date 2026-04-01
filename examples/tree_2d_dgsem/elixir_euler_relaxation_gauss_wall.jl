using OrdinaryDiffEqLowStorageRK
using OrdinaryDiffEq
using Trixi
using Plots

###############################################################################
# semidiscretization of the compressible Euler equations

eps = 1.0e-3
a = 1.0
equations = IncompressibleEulerRelaxationEquations2D(eps, a)

###############################################################################
# Get the DG approximation space
solver = DGSEM(polydeg = 4, surface_flux = flux_hll)

###############################################################################
# Mesh
coordinates_min = (-1.0, -1.0)
coordinates_max = (1.0, 1.0)

function initial_condition_gauss_wall(x, t, equations::IncompressibleEulerRelaxationEquations2D)
    p_eps = 0.0
    v1 = 2 * exp(-(x[1])^2 / 0.05)
    v2 = 0.0
    V_eps_11 = 0.0
    V_eps_12 = 0.0
    V_eps_21 = 0.0
    V_eps_22 = 0.0
    return SVector(p_eps, v1, v2, V_eps_11, V_eps_12, V_eps_21, V_eps_22)
end

mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 4,
                n_cells_max = 100_000,
                periodicity = false)

###############################################################################
# create the semi discretization object
boundary_conditions = (;
    x_neg = boundary_condition_slip_wall,
    x_pos = boundary_condition_slip_wall,
    y_neg = boundary_condition_do_nothing,
    y_pos = boundary_condition_do_nothing
)
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_gauss_wall, solver,
                                    source_terms = source_terms_constant,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 2.3*eps)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100
analysis_callback = AnalysisCallback(semi, interval = analysis_interval, analysis_integrals=())

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 10,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 0.3)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        #save_solution,
                        stepsize_callback
                        )

###############################################################################
# run the simulation
steps_vis = 20
#ssolver = KenCarp4(autodiff=false);
ssolver = CarpenterKennedy2N54(williamson_condition = false);
sol = solve(ode, ssolver;
            dt = tspan[2]/steps_vis/5, # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()..., callback = callbacks, saveat = range(ode.tspan..., length=steps_vis))
println("Simulation finished with code $(sol.retcode).")

doPlot = true
if !doPlot || sol.retcode != :Success
    println("Plotting skipped.")
else
    for i in 1:steps_vis
        pd = PlotData1D(sol.u[i], semi, slice=:x)
        p_v1 = plot(pd["v1"], ylim=(-2.0, 2.0))
        p_v2 = plot(pd["v2"], ylim=(-2.0, 2.0))
        p_peps = plot(pd["p_eps"])
        display(plot(p_v1, p_v2, p_peps, layout = (3, 1)))
    end
end

function plot_V_eps(sol, steps_vis)
    for i in 1:steps_vis
        pd = PlotData1D(sol.u[i], semi, slice=:x)
        p_veps_11 = plot(pd["V_eps_11"])
        p_veps_12 = plot(pd["V_eps_12"])
        p_veps_21 = plot(pd["V_eps_21"])
        p_veps_22 = plot(pd["V_eps_22"])
        display(plot(p_veps_11, p_veps_12, p_veps_21, p_veps_22, layout = (2, 2)))
    end
end