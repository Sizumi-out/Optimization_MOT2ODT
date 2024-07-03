# defenition of Holy Trait
export laser1064, fLaser1064, mkOptLat, ω_FORT
export MOT, fMOT
export constVal, givenFn

abstract type Lasers end
struct IsFn end
struct IsntFn end

abstract type CtrVals end
abstract type CtrFnStr <: CtrVals end

struct constVal <: CtrVals
    val::Real
end

struct givenFn <: CtrVals
    f::Function
end

include("Laser1064.jl")
include("MOT671.jl")

# Traits for fMOT, fLaser1064

IsFnStr(::fMOT) = IsFn()
IsFnStr(::fLaser1064) = IsFn()
IsFnStr(::fOptLat) = IsFn()
IsFnStr(::sOptLat) = IsFn()
IsFnStr(::Any) = IsntFn() # Default behavior is notFuncStr


IsFnStr(::Type{fLaser1064}) = IsFn() # IsFnStr(::Type{laser1064}) = IsntFn()

priority(::Dissipative) = Int8(1)
priority(::FORT) = Int8(2)

include("control_func.jl")  # struct for function of control value. include: FORT 

function get_val_t(str, sym, t = 0.0)
    get_val_t(IsFnStr(str), getfield(str, sym), t)
end
get_val_t(::IsFn, ValSource::givenFn, t) = ValSource.f(t)
get_val_t(::IsFn, ValSource::constVal, t) = ValSource.val
get_val_t(::IsFn, ValSource::T, t) where T <: CtrFnStr = δ_tras(t, ValSource)
get_val_t(::IsntFn, ValSource, t) = ValSource