export U_ODT, grad_U
export T0_Atom_onOL, T_dep_OL, T_from_ODT, ω_ODT_ax
export tθ0

"""tθ0(t, θ_ω, dt_int): struct for integral! called from functions in this file directly
* t::Float64
* θ_ω::Float64
* dt_int::Float64
"""
mutable struct tθ0
    t::Float64
    θ_ω::Float64
    dt_int::Float64
end

πcI(z, δ, P, lase) = P* (ω_FORT + δ)/lase.z_R/(1 + ((z - lase.z_f)/lase.z_R)^2)


## ODT func / not func
get_πcI(t, z, lase) = get_πcI(IsFnStr(lase), t, z, lase)
function get_πcI(::IsFn, t, z, lase)
    t_sft = t - lase.sect[1]
    return πcI(z, get_val_t(lase, :δ, t_sft), get_val_t(lase, :P, t_sft), lase) # return πcI(z, lase.δ(t_sft), lase.P(t_sft), lase)
end
get_πcI(::IsntFn, t, z, lase) = πcI(z, lase.δ, lase.P, lase)



# todo: when called paramsODT!, dt_int = 1e-6 cannot be set
"""
    paramsODT!(t, z, tθs, lase)

This function gets & updates the parameters of the ODT (Optical Dipole Trap) based on the given inputs.
## Arguments
- `t`: Time parameter
- `z`: Position parameter
- `tθs`: Time and angle parameter
- `lase`: Laser parameter

## Returns
πcI_in, πcI_ref, δ_in, 

"""
paramsODT!(t, z, tθs, lase) = paramsODT!(IsFnStr(lase), t, z, tθs, lase)

# dont integrate: because of not controling with function
function paramsODT!(::IsntFn, t, z, tθs, lase)
    δ_in  = lase.inc.δ
    δ_ref = lase.ref.δ
    πcI_in  = get_πcI(t, z, lase.inc, δ_in)
    πcI_ref = get_πcI(t, z, lase.ref, δ_ref)

    tθs.θ_ω = (δ_in - δ_ref)*t # θ_ω = Δω*t when ω doesnt dipend on time
    # println("t = $t in IsntFn: δ_in: $δ_in, \t δ_ref: $δ_ref")
    return πcI_in, πcI_ref, δ_in, δ_ref
end

# do integrate: because of controling with function
function paramsODT!(::IsFn, t, z, tθs, lase)
    t_Rshift = t - lase.sect[1]         # beginning of integral section 
    t_Lshift = tθs.t - lase.sect[1]     # terminal  of integral section

    δ_in, δ_ref = θ_by_ω!(tθs, lase, t_Lshift, t_Rshift)
    
    πcI_in  = get_πcI(t_Rshift, z, lase.inc, δ_in)
    πcI_ref = get_πcI(t_Rshift, z, lase.ref, δ_ref)
    # println("t = $t in IsFn: δ_in: $δ_in, \t δ_ref: $δ_ref")

    return πcI_in, πcI_ref, δ_in, δ_ref
end


"""
θ_by_ω!(tθs::tθ, lase, t_Lshift, t_Rshift)

Compute θ by ω.

Arguments
- `tθs::tθ`: The tθ object.
- `lase`: The laser object.
- `t_Lshift`: integral section of left side of the shifted time .
- `t_Rshift`: The right side of above.
"""
function θ_by_ω!(tθs::tθ, lase, t_Lshift, t_Rshift)
    δ_in  = get_val_t(lase.inc, :δ, t_Rshift) # δ_in  = get_δ(lase.inc, t_Rshift)
    δ_ref = get_val_t(lase.ref, :δ, t_Rshift) # δ_ref = get_δ(lase.ref, t_Rshift) 
    tθs.dt_int = t_Rshift - t_Lshift
    if tθs.dt_int > 0 # ⇔ if there is a itegral section
        # println("$t_Rshift : integrated dt $(tθs.dt_int)")
        t_m = (t_Lshift + t_Rshift)/2
        tθs.θ_ω += tθs.dt_int/6*(sympson(δ_in,  get_val_t(lase.inc, :δ, t_m), tθs.δ_L[1]) - sympson(δ_ref, get_val_t(lase.ref, :δ, t_m), tθs.δ_L[2]))
        # tθs.θ_ω += tθs.dt_int/2*((δ_in + tθs.δ_L[1]) - (δ_ref + tθs.δ_L[2]))        
        tθs.δ_L[1] = δ_in    # for next integral value
        tθs.δ_L[2] = δ_ref   # for next integral value
    # else (tθs.dt_int == 0) # ⇔ if there is no integral section => the integrated value is the same with before 
        # println("$t_Rshift : not integrated due to dt = $(tθs.dt_int)")
    end
    return δ_in, δ_ref
end

function θ_by_ω!(tθs::tθ0, lase, t_Lshift, t_Rshift)
    t_ax = collect(t_Lshift : tθs.dt_int : t_Rshift)
    δin_ary  = get_val_t.((lase.inc,), :δ, t_ax) # δin_ary  = get_δ.((lase.inc,), t_ax)
    δref_ary = get_val_t.((lase.ref,), :δ, t_ax) # δref_ary = get_δ.((lase.ref,), t_ax)
    tθs.θ_ω += (sympson(δin_ary) - sympson(δref_ary))*tθs.dt_int
    return δin_ary[end], δref_ary[end] # target is tθs.θ_ω 
end


sympson(leftVal, cVal, rVal) = leftVal + 4*cVal + rVal
sympson(ary::Vector{Float64}) = (ary[1] + ary[end] + 4*sum(ary[3:2:end-1]) + 2*sum(ary[2:2:end-1]))/3

## ODT
function U_ODT(t, z, lase1064::FORT, Atom, tθs0 = tθ0(0, 0, 0)) # tθs0 is not used 
    πcI = get_πcI(t, z, lase1064)
    # ω_L = ω_FORT + lase1064.δ
    # return -3*c/(2*Atom.ω_0^3)*Atom.γ*(1/(Atom.ω_0 - ω_L) + 1/(Atom.ω_0 + ω_L))*πcI
    return -Atom.U0*πcI
end # 

function ω_ODT_ax(lase1064::FORT, Atom)
    sqrt(2*abs(U_ODT(0, 0, lase1064, Atom))/ Atom.m )/lase1064.z_R
end

function T_from_ODT(lase1064, Atom)
    abs(U_ODT(0, 0, tθs0, lase1064, Atom))/k_b
end

## gradient
using ForwardDiff
function grad_Us(t, z, lase1064s, Atom, tθs)
    dU = 0.0
    for lase in lase1064s
        dU += ForwardDiff.derivative(z -> U_ODT(t, z, lase, Atom, tθs), z)[1]
    end
    return dU
end