export AryFORT


# fLaser δ_func dispatch vars: 5
function AryFORT(sects::Vector{Vector{N}} where N <: Real, δ_fn_ary::Vector{T} where T <: CtrFnStr, P, z_f::Real, z_R::Real)
    return [fLaser1064(sects[idx], δ_fn_ary[idx], P, z_f, z_R) for idx in eachindex(δ_fn_ary)]
end

# laser δ_func dispatch vars: 4
function AryFORT(δ_ary::Vector{T} where T <: CtrFnStr, P, z_f::Real, z_R::Real)
    return [laser1064(δ, P, z_f, z_R) for δ in δ_ary]
end

# fLaser z_R dispatch vars: 5
function AryFORT(sect, δ_func::CtrFnStr, P, z_f::Real, z_R_ary::Vector)
    return [fLaser1064(sect, δ_func, P, z_f, z_R) for z_R in z_R_ary]
end

# Laser z_R dispatch vars: 4
function AryFORT(δ, P, z_f::Real, z_R_ary::Vector)
    return [laser1064(δ, P, z_f, z_R) for z_R in z_R_ary]
end

# fLaser z_f dispatch vars: 5
function AryFORT(sect, δ_func::CtrFnStr, P, z_f_ary::Vector, z_R::Real)
    return [fLaser1064(sect, δ_func, P, z_f, z_R) for z_f in z_f_ary]
end

# Laser z_f dispatch vars: 4
function AryFORT(δ, P, z_f_ary::Vector, z_R::Real)
    return [laser1064(δ, P, z_f, z_R) for z_f in z_f_ary]
end