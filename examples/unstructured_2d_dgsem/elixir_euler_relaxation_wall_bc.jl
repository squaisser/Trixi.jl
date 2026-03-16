using OrdinaryDiffEqLowStorageRK
using OrdinaryDiffEq
using Trixi
using Plots

###############################################################################
eps = 1.0e-4
a = 1.0
equations = IncompressibleEulerRelaxationEquations2D(eps ,a)

@inline function uniform_flow_state(x, t, equations::IncompressibleEulerRelaxationEquations2D)

    # set the freestream flow parameters
    v_freestream = 0.3
    p_eps_freestream = 0.0

    theta = pi/90.0
    si, co = sincos(theta)
    v1 = v_freestream * co
    v2 = v_freestream * si

    V_eps_11 = 0.0
    V_eps_12 = 0.0
    V_eps_21 = 0.0
    V_eps_22 = 0.0

    return SVector(p_eps_freestream, v1, v2, V_eps_11, V_eps_12, V_eps_21, V_eps_22)
end

initial_condition = uniform_flow_state
source_terms = source_terms_constant

boundary_condition_uniform_flow = BoundaryConditionDirichlet(uniform_flow_state)
boundary_conditions = Dict(:Bottom => boundary_condition_uniform_flow,
                           :Top => boundary_condition_uniform_flow,
                           :Right => boundary_condition_uniform_flow,
                           :Left => boundary_condition_uniform_flow,
                           :Circle => boundary_condition_slip_wall)

###############################################################################
# Get the DG approximation space
solver = DGSEM(polydeg = 4, surface_flux = flux_hll)

###############################################################################
# Get the curved quad mesh from a file
mesh_file = Trixi.download("https://gist.githubusercontent.com/andrewwinters5000/8b9b11a1eedfa54b215c122c3d17b271/raw/0d2b5d98c87e67a6f384693a8b8e54b4c9fcbf3d/mesh_box_around_circle.mesh",
                           joinpath(@__DIR__, "mesh_box_around_circle.mesh"))

mesh = UnstructuredMesh2D(mesh_file)

###############################################################################
# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms,
                                    boundary_conditions = boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 5*eps)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval, analysis_integrals=())

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 10,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 0.03)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        #save_solution,
                        stepsize_callback
                        )

###############################################################################
# run the simulation

sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false);
            dt = 0.001,
            ode_default_options()..., callback = callbacks);
println("Simulation finished with code $(sol.retcode).")

doPlot = true
if !doPlot || sol.retcode != :Success
    println("Plotting skipped.")
else
    pd = PlotData2D(sol)
    p1 = plot(pd["v1"])
    plot!(getmesh(pd))
    display(p1)
    p2 = plot(pd["v2"])
    plot!(getmesh(pd))
    display(p2)
end