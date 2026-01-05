# The same setup as tree_1d_dgsem/elixir_euler_source_terms.jl
# to verify the StructuredMesh implementation against TreeMesh

using OrdinaryDiffEqLowStorageRK
using OrdinaryDiffEq
using Trixi
using Infiltrator

###############################################################################
# semidiscretization of the compressible Euler equations
#TODO: choose proper parameters
equations = IncompressibleEulerRelaxationEquations1D(0.001, 1.) #epsilon, a

initial_condition = initial_condition_riemann
# initial_condition = initial_condition_constant

boundary_condition_zero_dirichlet = BoundaryConditionDirichlet((x, t, equations) -> SVector(0.0, 0.0, 0.0))
boundary_conditions_hyperbolic = (;
                                  x_neg = boundary_condition_do_nothing, #BoundaryConditionDirichlet((x, t, equations) -> SVector(5.6e-3, 0.5, 2.5e-3)),
                                  x_pos = boundary_condition_do_nothing #BoundaryConditionDirichlet((x, t, equations) -> SVector(0.0, 1.0, 0.0))
                                 )

# Note that the expected EOC of 5 is not reached with this flux.
# Using `flux_hll` instead yields the expected EOC.

# Up to version 0.13.0, `max_abs_speed_naive` was used as the default wave speed estimate of
# `const flux_lax_friedrichs = FluxLaxFriedrichs(), i.e., `FluxLaxFriedrichs(max_abs_speed = max_abs_speed_naive)`.
# In the `StepsizeCallback`, though, the less diffusive `max_abs_speeds` is employed which is consistent with `max_abs_speed`.
# Thus, we exchanged in PR#2458 the default wave speed used in the LLF flux to `max_abs_speed`.
# To ensure that every example still runs we specify explicitly `FluxLaxFriedrichs(max_abs_speed_naive)`.
# We remark, however, that the now default `max_abs_speed` is in general recommended due to compliance with the 
# `StepsizeCallback` (CFL-Condition) and less diffusion.

solver = DGSEM(polydeg = 4, surface_flux = FluxLaxFriedrichs())

coordinates_min = (-1.0,)
coordinates_max = (1.0,)
cells_per_dimension = (16,)

#mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 5,
                n_cells_max = 30_000,
                periodicity = false
                )

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_homogeneous,
                                    boundary_conditions = boundary_conditions_hyperbolic
                                    )

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 0.0005)
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

stepsize_callback = StepsizeCallback(cfl = 0.8)

time_series = TimeSeriesCallback(semi, [(-0.5), (0.5)];
                                 interval=5,
                                 solution_variables=cons2cons,
                                 filename="tseries.h5")

callbacks = CallbackSet(summary_callback,   #works
                        analysis_callback,  #fixed: MethodError: no method matching cons2entropy if not analysis_integrals=() (default analysis integrals:entropy -> needs cons2entropy)
                        alive_callback,     #works
                        save_solution,      #works
                        #stepsize_callback,  #fixed: MethodError: no method matching max_abs_speeds
                        time_series         #works
                        )

###############################################################################
# run the simulation
sol = solve(ode, ImplicitEuler(autodiff=false);
            dt = 0.0001, # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()..., callback = callbacks);
println("Simulation finished with code $(sol.retcode).")

###############################################################################
# plot some results
using Plots
pd = PlotData1D(sol)

#analytical solution
p_eps_l, v1_l, V_eps_l = initial_condition(-0.5, 0.0, equations)
p_eps_r, v1_r, V_eps_r = initial_condition(0.5, 0.0, equations)

function analytical_solution_p_eps(x, t, equations::IncompressibleEulerRelaxationEquations1D, p_eps_l, v1_l, V_eps_l, p_eps_r, v1_r, V_eps_r)
    b = sqrt(equations.a + 1) / equations.epsilon
    if x/t < -b
        return p_eps_l
    elseif x/t < 0.0
        return p_eps_l + 1/(2*(equations.a + 1))*((p_eps_r-p_eps_l)+(V_eps_r-V_eps_l)) - 1/(2*b)*(v1_r - v1_l)
    elseif x/t < b
        return p_eps_r - 1/(2*(equations.a + 1))*((p_eps_r - p_eps_l)+(V_eps_r - V_eps_l)) - 1/(2*b)*(v1_r - v1_l)
    else
        return p_eps_r
    end
end
function analytical_solution_v1(x, t, equations::IncompressibleEulerRelaxationEquations1D, p_eps_l, v1_l, V_eps_l, p_eps_r, v1_r, V_eps_r)
    b = sqrt(equations.a + 1) / equations.epsilon
    if x/t < -b
        return v1_l
    elseif x/t < b
        return 1/2*(v1_r + v1_l) - b/(2*(equations.a + 1))*((p_eps_r - p_eps_l)+(V_eps_r - V_eps_l))
    else
        return v1_r
    end
end
function analytical_solution_V_eps(x, t, equations::IncompressibleEulerRelaxationEquations1D, p_eps_l, v1_l, V_eps_l, p_eps_r, v1_r, V_eps_r)
    b = sqrt(equations.a + 1) / equations.epsilon
    if x/t < -b
        return V_eps_l
    elseif x/t < 0.0
        return V_eps_l + equations.a/(2*(equations.a + 1))*((p_eps_r - p_eps_l)+(V_eps_r - V_eps_l)) - equations.a/(2*b)*(v1_r - v1_l)
    elseif x/t < b
        return V_eps_r - equations.a/(2*(equations.a + 1))*((p_eps_r - p_eps_l)+(V_eps_r - V_eps_l)) - equations.a/(2*b)*(v1_r - v1_l)
    else
        return V_eps_r
    end
end

p_eps_analytical = [analytical_solution_p_eps(x, sol.t[end], equations, p_eps_l, v1_l, V_eps_l, p_eps_r, v1_r, V_eps_r) for x in pd.x]
v1_analytical = [analytical_solution_v1(x, sol.t[end], equations, p_eps_l, v1_l, V_eps_l, p_eps_r, v1_r, V_eps_r) for x in pd.x]
V_eps_analytical = [analytical_solution_V_eps(x, sol.t[end], equations, p_eps_l, v1_l, V_eps_l, p_eps_r, v1_r, V_eps_r) for x in pd.x]
p1 = plot(pd["p_eps"], title = "Pressure Relaxation Variable", xlabel = "x", ylabel = "p_eps")#, ylims=(0.005, 0.006))
p1 = plot!(pd.x, p_eps_analytical, label = "analytical", lw=2, ls=:dash, color=:black)
p2 = plot(pd["v1"], title = "Velocity", xlabel = "x", ylabel = "v1")#, ylims=(0, 1))  
p2 = plot!(pd.x, v1_analytical, label = "analytical", lw=2, ls=:dash, color=:black)
p3 = plot(pd["V_eps"], title = "Relaxation Variable", xlabel = "x", ylabel = "V_eps")#, ylims=(0, 0.005))
p3 = plot!(pd.x, V_eps_analytical, label = "analytical", lw=2, ls=:dash, color=:black)
#display(plot(p1, p2, p3, layout = (3, 1)))
display(plot(p1))
display(plot(p2))
display(plot(p3))

#pd1 = PlotData1D(time_series, 1)
#pd2 = PlotData1D(time_series, 2)
#p1 = plot(pd1["p_eps"], label = "p_eps at x=-0.5", xlabel = "t", ylabel = "p_eps")
#plot!(pd2["p_eps"], label = "p_eps at x=0.5", xlabel = "t", ylabel = "p_eps")
#plot!(legend=:outerbottom, legendcolumns=2)
#display(plot(p1))
println("Plotting finished.")