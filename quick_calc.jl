export doppler_shift, zeeman_shift, T_D_Li, recoil_ene_Li
export get_w0, get_z_R
export T_from_v_Li, v_from_T_Li
export v_t_from_Δω, Δω_from_v_t
export LiAtom, NaAtom

LiAtom = select_atom("6Li")
NaAtom = select_atom("23Na")


λ = LiAtom.λ
mass = LiAtom.m

doppler_shift(v) =  v/λ

zeeman_shift(x, B) = μ_B*B*x

"Li: 速度を温度へ変換"
T_from_v_Li(v) = mass * v^2 /(2*k_b)

"Li: 温度を速度へ変換"
v_from_T_Li(T) = sqrt(2*k_b * T / mass)

"Get beam west from Rayleigh length"
get_w0(ω_L, z_R) = sqrt(2*c*z_R/ω_L)
get_w0(laser1064) = sqrt(2*c*laser1064.z_R/(laser1064.δ + ω_FORT))


"Get Rayleigh length from beam west"
get_z_R(ω_L, w0) = w0^2*ω_L/(2*c)

"ドップラー冷却限界"
T_D_Li = ħ*Li_γ/(2*k_b)

"recoil energy"
recoil_ene_Li = ħ^2*(2π/λ)^2 /(2*mass)

# useful funcs
T_D(Atom) = ħ*Atom.γ/(2*k_b)

ask_relaxetion_time_len(dt, Atom) = round(Int64, 1/Atom.γ/(2*pi*dt))
    
"速度を温度へ変換"
T_from_v(v, atom::Atom) = atom.m * v^2 / k_b
function T_from_v(v, atom::String) 
    Atom = select_atom(atom)
    return Atom.m * v^2 / k_b
end

"温度を速度へ変換"
v_from_T(T, Atom) = sqrt(k_b * T / Atom.m)