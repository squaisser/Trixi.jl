# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

struct IncompressibleEulerRelaxationEquations1D{RealT <: Real} <: 
        AbstractIncompressibleEulerRelaxationEquations{1, 3}
    epsilon::RealT  # relaxation parameter
    a::RealT        # viscosity parameter
end

function varnames(::typeof(cons2cons), ::IncompressibleEulerRelaxationEquations1D)
    return ("p_eps", "v1", "V_eps")
end

function initial_condition_constant(x, t, equations::IncompressibleEulerRelaxationEquations1D)
    #TODO: assign proper initial condition
    RealT = eltype(x)
    p_eps = convert(RealT, 1)
    v1 = convert(RealT, 1)
    V_eps = convert(RealT, 1)
    return SVector(p_eps, v1, V_eps)
end

function initial_condition_riemann(x, t, equations::IncompressibleEulerRelaxationEquations1D)
    #TODO: assign proper initial condition
    RealT = eltype(x)
    v1_l = convert(RealT, 0)
    v1_r = convert(RealT, 1)
    p_eps_l = convert(RealT, 1)*equations.epsilon^2
    p_eps_r = convert(RealT, 0.125)*equations.epsilon^2
    V_eps_l = convert(RealT, 1)*equations.epsilon^2
    V_eps_r = convert(RealT, 0.125)*equations.epsilon^2
    if x[1] < 0.0
        return SVector(p_eps_l, v1_l, V_eps_l)
    else
        return SVector(p_eps_r, v1_r, V_eps_r)
    end
end

@inline function source_terms_constant(u, x, t,
                              equations::IncompressibleEulerRelaxationEquations1D)
    p_eps, v1, V_eps = u
    s1 = zero(eltype(u))
    s2 = zero(eltype(u))
    s3 = v1^2 - (1/(equations.epsilon^2))*(V_eps)
    return SVector(s1, s2, s3)
end


@inline function flux(u, orientation::Integer,
                      equations::IncompressibleEulerRelaxationEquations1D)
    p_eps, v1, V_eps = u
    # Ignore orientation since it is always "1" in 1D
    f1 = v1
    f2 = 1/(equations.epsilon^2) * (V_eps + p_eps)
    f3 = equations.a*v1
    return SVector(f1, f2, f3)
end

@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                               equations::IncompressibleEulerRelaxationEquations1D)
    return sqrt(equations.a+1)/(equations.epsilon)
end

@inline function max_abs_speeds(u, equations::IncompressibleEulerRelaxationEquations1D)
    return (sqrt(equations.a+1)/(equations.epsilon),)
end

end # @muladd
