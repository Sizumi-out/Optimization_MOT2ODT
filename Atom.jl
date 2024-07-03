export Atom

include("TraitLasers.jl")   # only struct

const σplus = 1.0 # Left # What is the bst type? Float? Int8?
const σmns = -1.0 # Right # What is the bst type? Float? Int8?


""" The struct of Atoms
# Declare
Atom(m, λ, γ, γ_trns_streng_abs, γ_trns_streng_down)
* m: mass
* λ: wavelength
* γ: the natural width of transition [/rad]
* γ_trns_streng_abs: ある状態からの励起のときの遷移強度の比
* γ_trns_streng_down: ある状態からの光放出の遷移強度の比

# Additional fields
This struct has fields below in addiction to above. 
* ω_0: resonant frequency[/rad]
* k = 2π/λ: wavenumber
* U0: use for optical dipole trap_core
* v_k: recoil velocity
"""
struct Atom # const なので、型を指定する
    m::Float64
    λ::Float64
    ω_0::Float64
    k::Float64
    γ::Float64
    γ_trns_streng_abs::Vector{Tuple{Vector{Int16}, Vector{Int16}}} # Array(3) in Tuple(2) in Matrix(3)
    γ_trns_streng_down::Vector{Vector{Int16}} # Array(6) in Matrix(5)
    U0::Float64
    v_k::Float64
end


"原子の状態を指す構造体: term: ^Sなどの角軌道量, m_spin: 磁気量子数"
mutable struct AtomState # コンストラクタは要検証
    term::Rational{Int8} # S = 1, P = 2
    level::Rational{Int8} # 直接書く
    m_spin::Rational{Int8} # 直接書く
end


## contructor
function Atom(m, λ, γ, γ_trns_streng_abs, γ_trns_streng_down)
    ω_0 = ω_from_λ(λ)
    U0 = 3*c/(2*ω_0^3)*γ*(1/(ω_0 - ω_FORT) + 1/(ω_0 + ω_FORT))
    k = k_from_λ(λ)
    v_kick = ħ*k/m
    return Atom(m, λ, ω_from_λ(λ), k_from_λ(λ), γ, γ_trns_streng_abs, γ_trns_streng_down, U0, v_kick)
end


## State selecting functions
const sumγ = 270.0 # Transition Strengthの総和
const spin_under_down = [1//2, -1//2, 1//2, -1//2, 3//2] # 光放出時のスピンmの制御用配列
const level_under_down = [1//2, 1//2, 3//2, 3//2, 3//2] # 光放出時のFの制御用配列
const decide_trans = [1//4, 3//4, 9//4, 5//4, 15//4, 25//4] # どのγを用いるか制御する配列

abs_idx(state) = findfirst(isequal(state.level * abs(state.m_spin)), decide_trans)

abs_tpl_idx(state, pol) = Int((pol*sign(state.m_spin) + 3)/2) # σ+とスピンの情報よりTuple内(1or2)を指定

@views sum_strenth(Atom, γ_abs_idx, γ_abs_tpl) = sum(Atom.γ_trns_streng_abs[γ_abs_idx][γ_abs_tpl])

"光放出時の原子の量子状態の更新"
function down_atom_state!(AtomState, Atom)
    dicisison_state_val = rand()
    AtomState.term = 1//1
    γ_abs_idx = abs_idx(AtomState)
    distribution_func = 0.0
    cnt_dwn = UInt8(1)
    for dwn_streng in Atom.γ_trns_streng_down[γ_abs_idx]
        distribution_func += dwn_streng/sumγ # 累積分布
        if distribution_func <= dicisison_state_val 
            cnt_dwn += 1
        end
    end
    AtomState.m_spin = sign(AtomState.m_spin) * spin_under_down[6 - cnt_dwn]
    AtomState.level = level_under_down[6 - cnt_dwn]
    return # nothing
end

# todo: γ_abs_idx and γ_abs_tuple is multipliedly calculated with `laser-cooling`
"原子の量子状態の更新: (AtomState, Atom, polarization{1: σ+(left), -1: σ-(right)})"
function abs_atom_state!(AtomState, Atom, polarization, γ_abs_idx, γ_abs_tpl, norm_val)
    dicisison_state_val = rand()
    AtomState.term = 2//1
    # println("Reffered Array:", Atom.γ_trns_streng_abs[γ_abs_idx][γ_abs_tpl], ",\n  γ_abs_tuple:", γ_abs_tuple,",   level:", AtomState.level, ",  m_spin:", AtomState.m_spin)
    distribution_func = 0.0
    cnt = UInt8(0)
    for abs_streng in Atom.γ_trns_streng_abs[γ_abs_idx][γ_abs_tpl]
        distribution_func += abs_streng/norm_val # 累積分布
        if dicisison_state_val <= distribution_func
            cnt += 2 # Fは1/2なので、2足すと計算早いはず
        end
    end
    AtomState.level = (7 - cnt)//2
    AtomState.m_spin += polarization
    # println("   level:", AtomState.level, ",  m_spin:", AtomState.m_spin)
    return # nothing
end