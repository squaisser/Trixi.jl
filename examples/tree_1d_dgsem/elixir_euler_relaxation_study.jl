using OrdinaryDiffEqLowStorageRK
using OrdinaryDiffEq
using Trixi
using Infiltrator
using Plots
using LaTeXStrings

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

function plot_solution_riemann(sol, equations::IncompressibleEulerRelaxationEquations1D, combined=false)
    pd = PlotData1D(sol)
    
    #analytical solution
    p_eps_l, v1_l, V_eps_l = initial_condition_riemann(-0.5, 0.0, equations)
    p_eps_r, v1_r, V_eps_r = initial_condition_riemann(0.5, 0.0, equations)
    p_eps_analytical = [analytical_solution_p_eps(x, sol.t[end], equations, p_eps_l, v1_l, V_eps_l, p_eps_r, v1_r, V_eps_r) for x in pd.x]
    v1_analytical = [analytical_solution_v1(x, sol.t[end], equations, p_eps_l, v1_l, V_eps_l, p_eps_r, v1_r, V_eps_r) for x in pd.x]
    V_eps_analytical = [analytical_solution_V_eps(x, sol.t[end], equations, p_eps_l, v1_l, V_eps_l, p_eps_r, v1_r, V_eps_r) for x in pd.x]
    
    p1 = plot(pd["p_eps"], title = "Pressure Relaxation Variable", xlabel = "x", ylabel = "p_eps")#, ylims=(0.005, 0.006))
    plot!(pd.x, p_eps_analytical, label = "analytical", lw=2, ls=:dash, color=:black)
    p2 = plot(pd["v1"], title = "Velocity", xlabel = "x", ylabel = "v1")#, ylims=(0, 1))  
    plot!(pd.x, v1_analytical, label = "analytical", lw=2, ls=:dash, color=:black)
    p3 = plot(pd["V_eps"], title = "Relaxation Variable", xlabel = "x", ylabel = "V_eps")#, ylims=(0, 0.005))
    plot!(pd.x, V_eps_analytical, label = "analytical", lw=2, ls=:dash, color=:black)#, ylims=(-0.001, 0.004))
    
    if combined
        display(plot(p1, p2, p3, layout = (3, 1)))
    else
        display(plot(p1))
        display(plot(p2))
        display(plot(p3))
    end
end

function plot_solution_compare(sol_homogeneous, sol_nonhomogeneous, equations::IncompressibleEulerRelaxationEquations1D, combined=false)
    pdh = PlotData1D(sol_homogeneous)
    pdn = PlotData1D(sol_nonhomogeneous)
    
    #analytical solution
    p_eps_l, v1_l, V_eps_l = initial_condition_riemann(-0.5, 0.0, equations)
    p_eps_r, v1_r, V_eps_r = initial_condition_riemann(0.5, 0.0, equations)
    p_eps_analytical = [analytical_solution_p_eps(x, sol_homogeneous.t[end], equations, p_eps_l, v1_l, V_eps_l, p_eps_r, v1_r, V_eps_r) for x in pdh.x]
    v1_analytical = [analytical_solution_v1(x, sol_homogeneous.t[end], equations, p_eps_l, v1_l, V_eps_l, p_eps_r, v1_r, V_eps_r) for x in pdh.x]
    V_eps_analytical = [analytical_solution_V_eps(x, sol_homogeneous.t[end], equations, p_eps_l, v1_l, V_eps_l, p_eps_r, v1_r, V_eps_r) for x in pdh.x]
    
    p1 = plot(pdh.x, p_eps_analytical, label = "Analytical homogeneous", lw=2, ls=:dash, color=:black)
    plot!(pdh["p_eps"], xlabel = "x", label = "1D homogeneous")
    plot!(pdn["p_eps"], label = "1D nonhomogeneous", legend=:top, title="Pressure p_eps")
    

    p2 = plot(pdh.x, v1_analytical, label = "Analytical homogeneous", lw=2, ls=:dash, color=:black) 
    plot!(pdh["v1"], xlabel = "x", label = "1D homogeneous")
    plot!(pdn["v1"], label = "1D nonhomogeneous", legend=:topleft, title="Velocity u")
    
    
    p3 = plot(pdh.x, V_eps_analytical, label = "Analytical homogeneous", lw=2, ls=:dash, color=:black)
    plot!(pdh["V_eps"], xlabel = "x", label = "1D homogeneous")
    plot!(pdn["V_eps"], label = "1D nonhomogeneous", legend=:inside, title = "Relaxation variable V_eps")
    
    if combined
        display(plot(p1, p2, p3, layout = (3, 1)))
    else
        display(plot(p1))
        display(plot(p2))
        display(plot(p3))
    end
    return (p1, p2, p3)
end

###############################################################################
# semidiscretization of the compressible Euler equations
eps = 1.0e-3
a = 1.0
equations = IncompressibleEulerRelaxationEquations1D(eps, a)
tspan = (0.0, eps/sqrt(2)*0.5)

initial_condition = initial_condition_riemann

boundary_condition_zero_dirichlet = BoundaryConditionDirichlet((x, t, equations) -> SVector(0.0, 0.0, 0.0))
boundary_conditions_hyperbolic = (;
                                  x_neg = boundary_condition_do_nothing,
                                  x_pos = boundary_condition_do_nothing
                                 )

solver = DGSEM(polydeg = 4, surface_flux = FluxLaxFriedrichs())

coordinates_min = (-1.0,)
coordinates_max = (1.0,)

mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 5,
                n_cells_max = 30_000,
                periodicity = false
                )

semi_homogeneous = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_homogeneous,
                                    boundary_conditions = boundary_conditions_hyperbolic
                                    )
semi_nonhomogeneous = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    source_terms = source_terms_constant,
                                    boundary_conditions = boundary_conditions_hyperbolic
                                    )

###############################################################################
# ODE solvers, callbacks etc.

ode_homogeneous = semidiscretize(semi_homogeneous, tspan)
ode_nonhomogeneous = semidiscretize(semi_nonhomogeneous, tspan)

summary_callback = SummaryCallback()

analysis_interval = 100

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 100,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2cons)

stepsize_callback = StepsizeCallback(cfl = 0.3)

callbacks = CallbackSet(summary_callback,
                        alive_callback,
                        stepsize_callback
                        )

###############################################################################
# run the simulation
sol_homogeneous = solve(ode_homogeneous, CarpenterKennedy2N54(williamson_condition = false); dt = eps/10, ode_default_options()..., callback = callbacks);
println("Homogeneous simulation finished with code $(sol_homogeneous.retcode).")
sol_nonhomogeneous = solve(ode_nonhomogeneous, CarpenterKennedy2N54(williamson_condition = false); dt = eps/10, ode_default_options()..., callback = callbacks);
println("Nonhomogeneous simulation finished with code $(sol_nonhomogeneous.retcode).")

###############################################################################
# plot some results
p1, p2, p3 = plot_solution_compare(sol_homogeneous, sol_nonhomogeneous, equations, false)
path = "/home/ng143200/hiwi/Trixi.jl/plot/1D_explicit/"
savefig(p1, path*"p_eps.png")
savefig(p2, path*"v1.png")
savefig(p3, path*"V_eps.png")
println("Plotting finished.")