module AtomTrapCore

import Distributed.pmap

using ProgressMeter #, JLD2, DataFrames

export simul_kinetics, get_trap_final_DualAry, get_trap_final_MonoAry, eval_trap_prob, eval_trap_range, trap_dynamics, convey

include("constants.jl")
using .Constants
include("jld_io_write.jl")
using .JldWrite

include("trap.jl")
include("AlkalineAtom.jl")
include("quick_calc.jl")
include("struct_API.jl")
include("sampling.jl")

# sampling

""" trap_dynamics(dt, t_end, u0, Atom; sampling, v_length, direction)\n
# variables
## normal variables
* u = [x, v]: initial state of position and velocity \n
* dt: time step width,\t t_end: 終了時刻,\t 
* strAtom: define by  string like "6Li" or "23Na"\n
* TrapStructs: Varargs argument of struct Dissipative & FORT

## Keyword Arguments
* sampling: 配列のサンプリング間隔。デフォルトは1。配列長の約数であることが必須\n
rec_len: 速度の観測する飛ばしを勘案しない配列長。デフォルトでは時刻の配列数と同じ数。\n
            外から与えるときは平均を取る配列調にすれば良い。\n
            あるいは、round(Int64, (後ろから観測したい時間)/dt)\n
direction: apply potential mgz by true/false or giving to section to apply with `[[1e-3, 2e-3], [5e-3, 9e-3]]`
""" # ODT
function trap_dynamics(dt, t_end, u0, strAtom, TrapStructs...;
            sampling = 1, rec_len = t_end,
            gravity_sw = false, qs_rec = true, export_path = "N")

    Atom = select_atom(strAtom)

    # Progress always put valid.
    recorder = main(dt, t_end, u0, Atom, collect(TrapStructs), 
        sampling = sampling, rec_len = rec_len, 
        gravity_sw = gravity_sw, qs_rec_sw = qs_rec
    )
    export_path != "N" && jld_dynamics(recorder, export_path, dt, sampling, t_end, strAtom, TrapStructs, gravity_sw)
    # x_array = recorder.u_ary[1,:] # v_array = recorder.u_ary[2,:]
    if qs_rec
        return recorder.u_ary[1, :] , T_from_v.(recorder.u_ary[2, :], (Atom, )), recorder.u_ary[2, :], recorder.AS_ary
    else
        return recorder.u_ary[1, :], T_from_v.(recorder.u_ary[2, :], (Atom, )), recorder.u_ary[2, :]
    end
end


"return [x_initial, x_min, x_max, x_final, T_initial, T_min, T_max, T_final]" # ODT
function get_Mmx_v(dt, t_end, u0, Atom, TrapStructs...;
    sampling = 1, rec_len = t_end, gravity_sw = false)

    recorder = main(dt, t_end, u0, Atom, collect(TrapStructs), 
        sampling = sampling, rec_len = rec_len, 
       gravity_sw = gravity_sw, qs_rec_sw = false, prog_sw = false
   )
   x_array = first.(recorder.u_array)
   v_array = last.(recorder.u_array)
   arg_x_max = argmax(x_array)
   arg_x_min = argmin(x_array)

   return [u0[1], x_array[arg_x_min], x_array[arg_x_max], x_array[end], 
        T_from_v(u0[2], Atom), T_from_v(v_array[arg_x_min], Atom), T_from_v(v_array[arg_x_max], Atom), T_from_v(v_array[end], Atom)
    ]
end


""" get_trap_final(dt, t_end, u0, Atom, TrapStructs...; start_save = t_end, gravity_sw = false) \n
return [⟨|v|⟩, ⟨T⟩, |x[end]|]

# Arguments
* start_save_time: Time to start recording seconds before the end of the day.

"""
function get_trap_final(dt, t_end, u0, Atom, TrapStructs...; 
    start_save_time = t_end, gravity_sw = false)

    u_init = get_u0(u0, Atom)
    v_array_length = round(Int64, start_save_time/dt)
    recorder = main(dt, t_end, u_init, Atom, collect(TrapStructs), 
        rec_len = start_save_time, gravity_sw = gravity_sw, prog_sw = false, qs_rec_sw = false)
    # x_array = recorder.u_ary[1,:] # v_array = recorder.u_ary[2,:]
    return [u_init[1], u_init[2], sum(abs.(recorder.u_ary[2,:]))/v_array_length, sum(T_from_v.(recorder.u_ary[2,:], (Atom, )))/v_array_length, recorder.u_ary[1, end]]
end

get_u0(u0::Vector{Float64}, Atom) = u0
get_u0(boltzman::Boltzman, Atom) = Boltzmann_sampling(boltzman, Atom)



"""get_trap_final_ary(dt, t_end, u0_list, Atom, StructArray, TrapStructs...; 
        start_save_time = t_end, gravity_sw = false)
we can get [⟨v⟩, ⟨T⟩, |x[end]|] for each condition about u0_list & StructArray::Array{Lasers}
"""
function get_trap_final_DualAry(dt, t_end, u0_list, strAtom, StructArray, TrapStructs...;
        start_save_time = t_end, gravity_sw = false)
    
    Atom = select_atom(strAtom)
    results = []
    N = length(u0_list)*length(StructArray)
    println("Execute ArrayList mk_data_dep_detune per u0_list")
    p = Progress(N, dt=0.1,
            barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
            barlen=10)
    for u0 in u0_list
        push!(results, reduce(vcat, transpose.(
            pmap(laser -> get_trap_final(dt, t_end, u0, Atom, laser, TrapStructs..., 
            start_save_time = start_save_time, gravity_sw = gravity_sw), StructArray)))
        )
        next!(p)
    end
    ProgressMeter.finish!(p)
    return results
end

# Holy Traited Function
function get_trap_final_MonoAry(dt, t_end, u0, strAtom, SomeArray, TrapStructs...;
    start_save_time = t_end, gravity_sw=false, itrt=5)

    Atom = select_atom(strAtom)
    N = length(SomeArray)
    res_array = zeros(N, 3)
    println("Execute just δ dependency mk_data_dep_detune")
    p = Progress(N, dt=0.1,
                barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
                barlen=10)
    for _ in 1:itrt
        res_array += reduce(vcat, transpose.(
            get_trap_final_MonoAry(SomeArray[1], dt, t_end, u0, Atom, SomeArray, TrapStructs,
                start_save_time = start_save_time, gravity_sw = gravity_sw))
        )
        next!(p)
    end
    ProgressMeter.finish!(p)
    return res_array / itrt
end

# HT1: u0
function get_trap_final_MonoAry(::Array{Float64,1}, dt, t_end, u0, Atom, u0Array, TrapStructs;
    start_save_time = t_end, gravity_sw = false)
    return pmap(u0fA ->
            get_trap_final(dt, t_end, u0fA, Atom, TrapStructs, 
                start_save_time = start_save_time, gravity_sw = gravity_sw),
        u0Array)
end

# HT2: Lasers
function get_trap_final_MonoAry(::Lasers, dt, t_end, u0, Atom, LaserArray, TrapStructs;
    start_save_time = t_end, gravity_sw = false)
    return pmap(Laser ->
        get_trap_final(dt, t_end, u0, Atom, Laser, TrapStructs..., 
                start_save_time = start_save_time, gravity_sw = gravity_sw), 
        LaserArray)
end



"""`eval_trap_prob(dt, t_end, u0_list, Atom, StructArray, TrapStructs..., 
        start_save_time = t_end, gravity_sw = false, itrt = 50, threshold = 5.0e-3)`

2つのリストに対して、トラップ確率の評価を行う。同一条件で `itrt`回繰り返して判定する。\n
閾値`threshold`[m]で成功、失敗を判断
"""
function eval_trap_prob(dt, t_end, u0_list, strAtom, StructArray, TrapStructs...; 
        start_save_time = t_end, gravity_sw = false, itrt = 50, threshold = 5.0e-3)


    Atom = select_atom(strAtom)
    
    final_res_exp = [] # 離調の配列長分の平均 Vector[idx_δ] = Array[v, T, x]
    final_res_exp_pow2 = [] # 離調の配列長分の二乗平均 Vector[idx_δ] = Array[v, x].^2。
    prob_results = []# 離調の配列長分の平均 Vector[idx_δ] = Vector[idx_v0] = prob:: v0ごとの
    N = length(u0_list)*length(StructArray)
    p = Progress(N, dt=0.1,
        barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
        barlen=10)
    for laser in StructArray
        prob_array = zeros(length(u0_list))
        results = zeros(length(u0_list), 3) # 各v0でのほしい形（期待値格納）
        results_pow2 = zeros(length(u0_list), 3) # 各v0でのほしい形（二乗平均格納）
        for _ in 1:itrt # 初速度を変えて同一トラップで実施
            result = reduce(vcat, transpose.(
                   pmap(u0 -> get_trap_final(dt, t_end, u0, Atom, laser, TrapStructs..., 
                    start_save_time = start_save_time, gravity_sw = gravity_sw), 
                    u0_list)))
            # 初速度の配列長分のv, T, x つまり、pmap[idx_v] = [v, T, x]
            # reduce(vcat, transpose.(A)) # transposeで[a, b]→[a b]。これをvcatで縦に重ねて、reduceで集約
            results += result
            results_pow2 += result.^2
            prob_array += map(x -> abs(x) > threshold ? 0 : 1 , result[:, 3]) # 5mmでtrapを判定, 1と0でmap
            next!(p)
        end
        push!(final_res_exp, results./itrt) # δを変えた結果がpushされる。
        push!(final_res_exp_pow2, results_pow2./itrt) # δを変えたごとの二乗の平均の結果がpushされる。
        push!(prob_results, prob_array./itrt)
    end
    ProgressMeter.finish!(p)

    pow2_final_res_exp = similar(final_res_exp) # ⟨x⟩^2: 平均の二乗
    for idx = eachindex(final_res_exp)
        pow2_final_res_exp[idx] = final_res_exp[idx].^2
    end
    std_dev = similar(pow2_final_res_exp)
    for idx = eachindex(pow2_final_res_exp)
        std_dev[idx] = [sqrt.(dev) for dev in (final_res_exp_pow2[idx] - pow2_final_res_exp[idx])]
    end
    return final_res_exp, std_dev, prob_results
end


"""
return array of each (x, T) about [initial, min, max, final]
"""
function eval_trap_range(dt, t_end, u0_list, strAtom, trapStructs...;
    sampling = 1, gravity_sw = false)
    
    Atom = select_atom(strAtom)
    N = length(u0_list)
    println("Execute $(N) times simulation with the change of v0 at the same condition")
    p = Progress(N, dt=0.1,
            barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
            barlen=10)
    results_array = progress_pmap(u0 -> get_Mmx_v(dt, t_end, u0, Atom, trapStructs,
                        sampling = sampling, gravity_sw = gravity_sw), 
                        u0_list, progress=p)
    return results_array
end


function simul_kinetics(dt, t_end, u0, strAtom, SomeArray, TrapStructs...;
        width = length(SomeArray), shift = 0, sampling_rate = 1, gravity_sw = false, qs_rec = true)
    
    N = length(SomeArray)
    p = Progress(N, dt=10,
        barglyphs=BarGlyphs('|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',),
        barlen=10)
    results = simul_kinetics(SomeArray[1], p, dt, t_end, u0, strAtom, SomeArray, TrapStructs, 
        width = width, shift = shift, sampling_rate = sampling_rate, gravity_sw = gravity_sw, qs_rec = qs_rec
    )
    return results
end

function simul_kinetics(::Array{Float64,1}, p, dt, t_end, u0, strAtom, u0Array, TrapStructs;
        width = length(u0Array), shift = 0, sampling_rate = 1, gravity_sw = false, qs_rec = true)
    return progress_pmap(u0 -> trap_dynamics(dt, t_end, u0, Atom, TrapStructs...,
            sampling = sampling_rate, gravity_sw = gravity_sw, qs_rec = qs_rec), 
        u0Array[1 + shift:width + shift], progress=p # 設定値が変更される
    )
end

function simul_kinetics(::T, p, dt, t_end, u0, strAtom, laserArray, TrapStructs...;
        width = length(laserArray), shift = 0, sampling_rate = 1, gravity_sw = false, qs_rec = true)  where T <: Lasers
    return progress_pmap(laser -> trap_dynamics(dt, t_end, u0, Atom, laser, TrapStructs...,
            sampling = sampling_rate, gravity_sw = gravity_sw, qs_rec = qs_rec), 
        laserArray[1 + shift:width + shift], progress=p # 設定値が変更される
    )
end

function get_terminal_time(structs, t_end)
    filter!(str->IsFnStr(str) == IsFn(), structs)
    return maximum(vcat([str.sect[2] for str in structs], t_end)) 
end

end # end module