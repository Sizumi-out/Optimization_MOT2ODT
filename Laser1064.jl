abstract type FORT <: Lasers end # if these lasers are ovewrapped, these make stading wave

const λ_FORT = 1064e-9          # 1064nm
const ω_FORT = ω_from_λ(λ_FORT) # 1064nmの角周波数
const k_FORT = k_from_λ(λ_FORT) # 1064nmの波数


"""laser1064(δ, P, z_f, z_R): struct for 1064nm trap: 
# variables
* δ::Float64    : detuned frequency by AOM from 1064nm
* P::Float64    : power
* z_f::Float64  : focus coodinate
* z_R::Float64  : Rayleigh range
# how to use
* use with `const ω_FORT = 2*pi*c/1064e-9`
"""
struct laser1064 <:FORT
    δ::Float64 # AOMによる離調
    P::Float64 # power
    z_f::Float64 # focus coodinate
    z_R::Float64 # Rayleigh range
end


"""1064nm laserの構造体: laser1064(δ, P, z_f, z_R)
# variables
* sect::Vector{Float64}: give applied section this is used. example `[20, 21]`
* δ  <: CtrVals    : detuned frequency by AOM from 1064nm
* P  <: CtrVals    : power
* z_f:: Float64  : focus coodinate
* z_R:: Float64  : Rayleigh range

# how to use
* example
    * funcs : `10exp(-5*t)`
    * sect  : `[1e-3, 4e-3]`
* use with `const ω_FORT = 2*pi*c/1064e-9`
"""
struct fLaser1064{T <: CtrVals, S <: CtrVals} <:FORT
    sect::Vector{Float64}
    δ::T # AOMによる離調
    P::S # power
    z_f::Float64 # focus coodinate
    z_R::Float64 # Rayleigh range
end

function fLaser1064(sect, δ_put, P_put, z_f, z_R)
    δ = valJudge(δ_put)
    P = valJudge(P_put)
    return fLaser1064(sect, δ, P, z_f, z_R)
end

valJudge(val::Real) = constVal(val)
valJudge(val::Function) = givenFn(val)
valJudge(val::CtrFnStr) = val