# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

struct IncompressibleEulerRelaxationEquations2D{RealT <: Real} <: 
        AbstractIncompressibleEulerRelaxationEquations{2, 7}
    epsilon::RealT  # relaxation parameter
    a::RealT        # viscosity parameter
end

function varnames(::typeof(cons2cons), ::IncompressibleEulerRelaxationEquations2D)
    return ("p_eps", "v1", "v2", "V_eps_11", "V_eps_12", "V_eps_21", "V_eps_22")
end

function initial_condition_constant(x, t,
                                  equations::IncompressibleEulerRelaxationEquations2D)
    #TODO: assign proper initial condition
    RealT = eltype(x)
    p_eps = convert(RealT, 1)
    v1 = convert(RealT, 1)
    v2 = convert(RealT, 0)
    V_eps_11 = convert(RealT, 1)
    V_eps_12 = convert(RealT, 0)
    V_eps_21 = convert(RealT, 0)
    V_eps_22 = convert(RealT, 0)
    return SVector(p_eps, v1, v2, V_eps_11, V_eps_12, V_eps_21, V_eps_22)
end

function initial_condition_riemann(x, t,
                                  equations::IncompressibleEulerRelaxationEquations2D)
    #TODO: assign proper initial condition
    RealT = eltype(x)
    v1_l = convert(RealT, 0)
    v1_r = convert(RealT, 1)
    v2_l = convert(RealT, 0)
    v2_r = convert(RealT, 0)
    p_eps_l = convert(RealT, 1)*equations.epsilon^2
    p_eps_r = convert(RealT, 0.125)*equations.epsilon^2
    V_eps_11_l = convert(RealT, 1)*equations.epsilon^2
    V_eps_11_r = convert(RealT, 0.125)*equations.epsilon^2
    V_eps_12_l = convert(RealT, 0)
    V_eps_12_r = convert(RealT, 0)
    V_eps_21_l = convert(RealT, 0)
    V_eps_21_r = convert(RealT, 0)
    V_eps_22_l = convert(RealT, 0)*equations.epsilon^2
    V_eps_22_r = convert(RealT, 0)*equations.epsilon^2
    if x[1] < 0.0
        return SVector(p_eps_l, v1_l, v2_l, V_eps_11_l, V_eps_12_l, V_eps_21_l, V_eps_22_l)
    else
        return SVector(p_eps_r, v1_r, v2_r, V_eps_11_r, V_eps_12_r, V_eps_21_r, V_eps_22_r)
    end
end

@inline function source_terms_homogeneous(u, x, t,
                              equations::IncompressibleEulerRelaxationEquations2D)
    s = zero(eltype(u))
    return SVector(s, s, s, s, s, s, s)
end

@inline function source_terms_constant(u, x, t,
                              equations::IncompressibleEulerRelaxationEquations2D)
    p_eps, v1, v2, V_eps_11, V_eps_12, V_eps_21, V_eps_22 = u
    s1 = zero(eltype(u))
    s2 = zero(eltype(u))
    s3 = zero(eltype(u))
    s4 = v1^2 - 1/(equations.epsilon^2) * V_eps_11
    s5 = v1*v2 - 1/(equations.epsilon^2) * V_eps_12
    s6 = v2*v1 - 1/(equations.epsilon^2) * V_eps_21
    s7 = v2^2 - 1/(equations.epsilon^2) * V_eps_22
    return SVector(s1, s2, s3, s4, s5, s6, s7)
end

@inline function flux(u, orientation::Integer,
                      equations::IncompressibleEulerRelaxationEquations2D)
    p_eps, v1, v2, V_eps_11, V_eps_12, V_eps_21, V_eps_22 = u
    if orientation == 1
        f1 = v1
        f2 = 1/(equations.epsilon^2) * (V_eps_11 + p_eps)
        f3 = 1/(equations.epsilon^2) * V_eps_12
        f4 = equations.a * v1
        f5 = 0
        f6 = equations.a * v2
        f7 = 0
    elseif orientation == 2
        f1 = v2
        f2 = 1/(equations.epsilon^2) * V_eps_21
        f3 = 1/(equations.epsilon^2) * (V_eps_22 + p_eps)
        f4 = 0
        f5 = equations.a * v1
        f6 = 0
        f7 = equations.a * v2
    else
        throw("Invalid orientation: $orientation")
    end
    return SVector(f1, f2, f3, f4, f5, f6, f7)
end

@inline function flux(u, normal_direction::AbstractVector,
                    equations::IncompressibleEulerRelaxationEquations2D)
    p_eps, v1, v2, V_eps_11, V_eps_12, V_eps_21, V_eps_22 = u
    nx, ny = normal_direction
    v_normal = v1*nx + v2*ny

    f1 = v_normal
    f2 = 1/(equations.epsilon^2) * (V_eps_11*nx + V_eps_21*ny + p_eps*nx)
    f3 = 1/(equations.epsilon^2) * (V_eps_12*nx + V_eps_22*ny + p_eps*ny)
    f4 = equations.a * v1 * nx
    f5 = equations.a * v1 * ny
    f6 = equations.a * v2 * nx
    f7 = equations.a * v2 * ny
    return SVector(f1, f2, f3, f4, f5, f6, f7)
end

@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                               equations::IncompressibleEulerRelaxationEquations2D)
    return sqrt(equations.a+1)/(equations.epsilon)
end

@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                               equations::IncompressibleEulerRelaxationEquations2D)
    return sqrt(equations.a+1)/(equations.epsilon)
end

@inline function max_abs_speeds(u, equations::IncompressibleEulerRelaxationEquations2D)
    return (sqrt(equations.a+1)/(equations.epsilon), sqrt(equations.a+1)/(equations.epsilon))
end

end # @muladd