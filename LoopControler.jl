const dim2 = Int8(2)

struct fCtrler
    swMat::Matrix{Bool}
    t_ax_len::Int64
    not_rec_len::Int64
    sampling::Int64
    prog_sw::Bool
end

struct Ctrler
    t_ax_len::Int64
    not_rec_len::Int64
    sampling::Int64
    prog_sw::Bool
end

mutable struct Recorder
    u_ary::Matrix{Float64}
    AS_ary::Matrix{fieldtype(AtomState, 1)}
    cnt::Int64
end

mutable struct RecorderGo
    u_ary::Matrix{Float64}
    cnt::Int64
end

mutable struct tvx
    u::MVector{2, Float64} # todo: MVector → Vector
    cnt::Int64
end

"""tgFn: time and gravity
* dt::Float64
* dt_half::Float64
* acc_ary::Vector{Bool}
"""
struct tgFn
    dt::Float64
    dt_half::Float64
    g_sw_ary::Vector{Bool}
end

"""tgFn: time and gravity
* dt::Float64
* dt_half::Float64
* g::Float64
"""
struct tg
    dt::Float64
    dt_half::Float64
    g::Float64
end

"""struct for integral! called from APIs (acc is only irrelevant value)
"""
mutable struct tθ
    t::Float64
    θ_ω::Float64
    δ_L::Vector{Float64}
    dt_int::Float64
    acc::Float64
end


## function for Loop structs (constructer)
"""mk_Recorder_array(field_num, N)
* field_num:: 1 => u_ary, 2 => AS_ary
* N:: length of array
"""
mk_Recorder_array(field_num::Integer, N) = mk_Recorder_array(eltype(fieldtype(Recorder, field_num)), N)
mk_Recorder_array(::Type{Float64}, N) = zeros(Float64, dim2, N)
mk_Recorder_array(::Type{Rational{Int8}}, N) = zeros(Rational{Int8}, fieldcount(AtomState), N)

"decide gravity(true/false): make gravity array"
gravity(grav_sw) = grav_sw ? gra_acc : 0.0

# for optical lattice

# external constructors
function mktθ(trapStructs, tg)  # need to fix: current = apply fucntion from dt ~ t_end 
    get_fLat = map(trapStr ->typeof(trapStr) == fOptLat, trapStructs)
    if any(get_fLat)
        println("set initial tθ")
        fLat = trapStructs[get_fLat][1] # decide Δω from 1st fOptlat
        δinit = [get_val_t(fLat.inc, :δ, tg.dt), get_val_t(fLat.ref, :δ, tg.dt)]
        return tθ(0.0, 0.0, δinit, 0.0, 0.0)
    else
        return tθ(0.0, 0.0, [0.0, 0.0], 0.0, 0.0) # not use [0.0, 0.0]
    end
end

function Recorder(v_array_len)
    u_array = mk_Recorder_array(1, v_array_len) # x, v
    atom_state_array = mk_Recorder_array(2, v_array_len)
    return Recorder(u_array, atom_state_array, 1)
end

function RecorderGo(v_array_len)
    u_array = mk_Recorder_array(1, v_array_len) # x, v
    return RecorderGo(u_array, 1)
end

tgFn(dt, g_sw) = tgFn(dt, dt/2, g_sw)
tg(dt, g_sw) = tg(dt, dt/2, g_sw)


# struct generater
function geneRecorder(v_array_len, qs_rec_sw)
    if qs_rec_sw
        return Recorder(v_array_len)
    else
        return RecorderGo(v_array_len)
    end
end

function ErrLoopStructs(trapStructs, t_ax_len, t_rec_len)
    judge = IsFnStr.(trapStructs)
    if !allequal(judge)
        throw(DomainError("Struct with different control logic are given as variables\n Structs: $judge"))
    elseif t_ax_len < t_rec_len # 速度を保存する時間の長さ
        throw(DomainError("t_rec_len should be shorter than t_axis_len"))
    end
end

""" LoopStructs(dt, t_end, sampling, v_length, g_sects, prog_sw, Structs...)
Dual constructor for struct `Ctrler` & `Recorder`
""" 
function LoopStructs(dt, t_end, sampling, t_rec_len, g_sects, qs_rec_sw, prog_sw, trapStructs) # Trait
    t_axis_len = 1 + round(Int64, t_end/dt) # 時間軸ベクトルの長さ
    v_length = 1 + round(Int64, t_rec_len/dt)
    v_array_len = ceil(Int64, v_length/sampling) # 速度の配列長

    ErrLoopStructs(trapStructs, t_axis_len, t_rec_len)

    not_record_len = t_axis_len - v_length
    # println("t_axis_len: $t_axis_len,  v_array_len: $v_array_len,  v_length: $v_length, not_record_len: $not_record_len")
    t_axis = collect(range(0, t_end, length=t_axis_len))
    # println("t_axis_len: $t_axis_len,  t_axis end: $(t_axis[end])")
    tgStr = decide_grav(t_axis, dt, g_sects, t_axis_len)
    LoopStructs(IsFnStr(trapStructs[1]), sampling, tgStr, qs_rec_sw, prog_sw, trapStructs, 
        t_axis, t_axis_len, v_array_len, not_record_len)
end

function LoopStructs(::IsFn, sampling, tgStr, qs_rec_sw, prog_sw, trapStructs, 
                    t_axis, t_axis_len, v_array_len, not_record_len)
    swMat = mk_swMat(trapStructs, t_axis, t_axis_len)
    return fCtrler(swMat, t_axis_len, not_record_len, sampling, prog_sw), geneRecorder(v_array_len, qs_rec_sw), tgStr
end

function LoopStructs(::IsntFn, sampling, tgStr, qs_rec_sw, prog_sw, trapStructs, # no one knows what struct is trapStructs
                    t_axis, t_axis_len, v_array_len, not_record_len)
    return Ctrler(t_axis_len, not_record_len, sampling, prog_sw), geneRecorder(v_array_len, qs_rec_sw), tgStr
end


function decide_grav(t_axis, dt, g_sects::Any, t_axis_len) 
    g_sw_ary = mk_sw_ary(t_axis, g_sects, t_axis_len)
    return tgFn(dt, dt/2, g_sw_ary)
end
function decide_grav(t_axis, dt, g_sw::Bool, t_axis_len) # g_sects = true/false
    g = gravity(g_sw)
    return tg(dt, dt/2, g)
end



## Make switching matrix
"""
mk_sw_ary(t_axis, sect...): make a true/false array
* t_axis::  Vector{Float64}
* sect: input section to turn on gravity
    * Vector{Float64}, Vector{Float64}, Vector{Float64} ⋯⋯
"""
function mk_sw_ary(t_axis, sect::Vector{T}, t_axis_len) where T <: Number
    sw_vec = fill(false, t_axis_len)
    sw_vec = map(t-> sect[1] <= t <= sect[2], t_axis)
    return sw_vec
end

function mk_sw_ary(t_axis, sects::Vector{Vector{T}}, t_axis_len) where T <: Number
    sw_vec = fill(false, t_axis_len)
    sw_vec = map(t-> any(isInSect.(sects, t)), t_axis)
    return sw_vec
end
isInSect(sect, t) = sect[1] <= t <= sect[2]


"""
mk_switches(fstructs, t_axis, sect): make a switching array to turn on or off ODT MOT
"""
function mk_swMat(fstructs, t_axis, t_axis_len)
    sw_Mat = Matrix{Bool}(undef, length(fstructs), t_axis_len)
    for (i, fstr) in enumerate(fstructs)
        sw_Mat[i, :] = mk_sw_ary(t_axis, fstr.sect, t_axis_len)
    end
    return sw_Mat
end


## functions
# todo: check what is fast
function record!(recorder::Recorder, tvx, atomstate)
    @views recorder.u_ary[:, recorder.cnt] = tvx.u
    recorder.AS_ary[1, recorder.cnt] = atomstate.term
    recorder.AS_ary[2, recorder.cnt] = atomstate.level
    recorder.AS_ary[3, recorder.cnt] = atomstate.m_spin
    recorder.cnt += 1
    return # nothing
end

function record!(recorder::RecorderGo, tvx, atomstate)
    @views recorder.u_ary[:, recorder.cnt] = tvx.u
    recorder.cnt += 1
    return # nothing
end


function judge_record(tvx, ctrler)
    len_to_record =  tvx.cnt - ctrler.not_rec_len
    return (len_to_record > 0) && ((len_to_record - 1) % ctrler.sampling == 0) # `-1` is put because the first index of array is 1 
end