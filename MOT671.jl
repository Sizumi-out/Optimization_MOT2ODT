abstract type Dissipative <: Lasers end

"""MOT(δ, s, B): struct for MOT 
# variables
* δ::Float64    : detuned frequency from ⁶Li. this should be inputed more than 0 because positive val is deeled as red detune.
* s::Float64    : = I/I_sat strength of laser
* B::Float64    : magnetic field gradient
"""
struct MOT <: Dissipative
    δ::Float64 # 共鳴からの離調
    s::Float64 # = I/I_sat: strength of laser
    B::Float64 # magnetic field gradient
end


"""fMOT(δ, δSect, s, sSect, B, Bsect): struct for MOT to control functionally
# variables
* sect::Vector{Float64} : give applied fMOT section. example `[0, 20]`
* δ::Function    : detuned frequency from ⁶Li. this should be inputed more than 0 because positive val is deeled as red detune.
* s::Function    : = I/I_sat strength of laser
* B::Function    : magnetic field gradient

# example
* funcs : `10exp(-5*t)`
* sect  : `[1e-3, 4e-3]`
"""
struct fMOT{T <: CtrVals} <: Dissipative
    sect::Vector{Float64}
    δ::T # 共鳴からの離調
    s::T # = I/I_sat: strength of laser
    B::T # magnetic field gradient
end