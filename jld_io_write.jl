"""
Export by Jld(binary) File. \n\n
If you want export data, excute following code.
    
    using DataFrames, JLD \n
    data = load(path*"01single_dynamics.jld2")
    
"""
module JldWrite

# CSV を吐く関数も作りたい
using JLD2
using Base.Threads

export jld_dynamics, jld_convey
# export jld_write_laser_cooling, jld_write_mk_single_array, jld_write_mk_data_dep_detune, jld_write_eval_trap_prob

function jld_dynamics(recorder, path, dt, sampl, t_end, strAtom, TrapStructs, gravity_sw)
    println("Export your input and the result by jld file")
    save(path*"/dynamics_$(strAtom).jld2", "Result", recorder,
        "t_ax", collect(0:sampl*dt:t_end), "dt", dt, 
        "TrapStructs", TrapStructs, "gravity_sw", gravity_sw
    )
    return #nothing
end

function jld_convey(path, res, strSome, SomeArray, t_end, boltzman, strAtom, TrapStructs...;
    N = 500, gravity_sw = false)

    save_path = path*"final_val_$(strAtom)_$(N)times.jld2"
    println("Export your input and the result by jld file as $save_path")
    save(save_path, "Result", res, "$(strSome)", SomeArray, 
        "t_end", t_end, "boltzman", boltzman,
        "TrapStructs", TrapStructs, "gravity_sw", gravity_sw
    )
    return #nothing
end
# widthとshiftを駆使して、ファイル名を設定
# δのディスパッチ: function
# function jld_write_mk_single_array(results, dt, t_end, s::Real, δ_list::Array{Function, 1}, 
#     B::Real, ω_p, width, shift; sampling_rate = 1, direction = "none"::String, path = pwd())

#     println("Export your input and the result by jld file which manipulated δ function")

#     t_axis_step = collect(0:sampling_rate *dt:t_end )
#     df_array= Array{DataFrame}(undef, width)
#     @threads for idx = eachindex(results) # output
#         df_array[idx] = hcat(DataFrame(time=t_axis_step, v=results[idx][1], 
#             T = results[idx][2], x=results[idx][3],), 
#             DataFrame(results[idx][4], ["term", "level", "state"])
#         )
#     end
#     filename = "/02mult_dynamics_δ_function.jld2"
#     save(path*filename, 
#         "df_Result", df_array,
#         "δ_list", δ_list, 
#         "dt", dt,   "t_end", t_end, 
#         "s", s,     "B", B,     "ω_p", ω_p,
#         "sampling_rate", sampling_rate,
#         "direction", direction, 
#         "width", width, "shift", shift
#     )
#     return # nothing
# end


# # δのディスパッチ
# function jld_write_mk_single_array(results, dt, t_end, s::Real, δ_list::Array{Float64, 1}, 
#     B::Real, ω_p, width, shift; sampling_rate = 1, direction = "none"::String, path = pwd())

#     println("Export your input and the result by jld file")

#     t_axis_step = collect(0:sampling_rate *dt:t_end )
#     df_array= Array{DataFrame}(undef, width)
#     @threads for idx = eachindex(results) # output
#         df_array[idx] = hcat(DataFrame(time=t_axis_step, v=results[idx][1], 
#             T = results[idx][2], x=results[idx][3],), 
#             DataFrame(results[idx][4], ["term", "level", "state"])
#         )
#     end
#     bgn_str   =    string(round(δ_list[begin + shift]/(2π*1e6); sigdigits = 2))
#     terminal_str = string(round(δ_list[width + shift]/(2π*1e6); sigdigits = 2))
#     filename = "/02mult_dynamics_δ"*bgn_str*"-"*terminal_str*"MHz.jld2"
#     save(path*filename, 
#         "df_Result", df_array,
#         "δ_list", δ_list, 
#         "dt", dt,   "t_end", t_end, 
#         "s", s, "B", B, "ω_p", ω_p,
#         "sampling_rate", sampling_rate,
#         "direction", direction, 
#         "width", width, "shift", shift
#     )
#     return # nothing
# end


# # s のディスパッチ
# function jld_write_mk_single_array(results, dt, t_end, s_list::Array{Float64, 1}, δ::Real, 
#     B::Real, ω_p, width, shift; sampling_rate = 1, direction = "none"::String, path = pwd())

#     println("Export your input and the result by jld file")

#     t_axis_step = collect(0:sampling_rate *dt:t_end )
#     df_array= Array{DataFrame}(undef, width)
#     @threads for idx = eachindex(results) # output
#         df_array[idx] = hcat(DataFrame(time=t_axis_step, v=results[idx][1], 
#             T = results[idx][2], x=results[idx][3],), 
#             DataFrame(results[idx][4], ["term", "level", "state"])
#         )
#     end
#     bgn_str   =    string(round(s_list[begin + shift]; sigdigits = 2))
#     terminal_str = string(round(s_list[width + shift]; sigdigits = 2))
#     filename = "/03mult_dynamics_s_"*bgn_str*"-"*terminal_str*".jld2"
#     save(path*filename, 
#         "df_Result", df_array,
#         "s_list", s_list, 
#         "dt", dt,   "t_end", t_end, 
#            "δ", δ, "B", B, "ω_p", ω_p,
#         "sampling_rate", sampling_rate,
#         "direction", direction, 
#         "width", width, "shift", shift
#     )
#     return # nothing
# end
# # s のディスパッチ: function
# function jld_write_mk_single_array(results, dt, t_end, v0, s_list::Array{Real, 1}, δ::Function, 
#     B::Real, ω_p, width, shift; sampling_rate = 1, direction = "none"::String, path = pwd())

#     println("Export your input and the result by jld file which manipulated δ function")

#     t_axis_step = collect(0:sampling_rate *dt:t_end )
#     df_array= Array{DataFrame}(undef, width)
#     @threads for idx = eachindex(results) # output
#         df_array[idx] = hcat(DataFrame(time=t_axis_step, v=results[idx][1], 
#             T = results[idx][2], x=results[idx][3],), 
#             DataFrame(results[idx][4], ["term", "level", "state"])
#         )
#     end
#     filename = "/03mult_dynamics_s_δfunction.jld2"
#     save(path*filename, 
#         "df_Result", df_array,
#         "s_list", s_list, 
#         "dt", dt,   "t_end", t_end, 
#         "v0", v0,   "δ", δ, "B", B, "ω_p", ω_p,
#         "sampling_rate", sampling_rate,
#         "direction", direction, 
#         "width", width, "shift", shift
#     )
#     return # nothing
# end



# # B のディスパッチ
# function jld_write_mk_single_array(results, dt, t_end, s::Real, δ::Real, B_list::Array{Float64, 1}, 
#     ω_p, width, shift; sampling_rate = 1, direction = "none"::String, path = pwd())

#     println("Export your input and the result by jld file")

#     t_axis_step = collect(0:sampling_rate *dt:t_end )
#     df_array= Array{DataFrame}(undef, width)
#     @threads for idx = eachindex(results) # output
#         df_array[idx] = hcat(DataFrame(time=t_axis_step, v=results[idx][1], 
#             T = results[idx][2], x=results[idx][3],), 
#             DataFrame(results[idx][4], ["term", "level", "state"])
#         )
#     end
#     bgn_str   =    string(B_list[begin + shift])
#     terminal_str = string(B_list[width + shift])
#     filename = "/04mult_dynamics_B_"*bgn_str*"-"*terminal_str*".jld2"
#     save(path*filename, 
#         "df_Result", df_array,
#         "B_list", B_list,
#         "dt", dt,   "t_end", t_end, 
#            "δ", δ, "s", s, "ω_p", ω_p,
#         "sampling_rate", sampling_rate,
#         "direction", direction, 
#         "width", width, "shift", shift
#     )
#     return # nothing
# end
# # B のディスパッチ: Function
# function jld_write_mk_single_array(results, dt, t_end, v0, s::Real, δ::Function, B_list::Array{Float64, 1}, 
#     ω_p, width, shift; sampling_rate = 1, direction = "none"::String, path = pwd())

#     println("Export your input and the result by jld file which manipulated δ function")

#     t_axis_step = collect(0:sampling_rate *dt:t_end )
#     df_array= Array{DataFrame}(undef, width)
#     @threads for idx = eachindex(results) # output
#         df_array[idx] = hcat(DataFrame(time=t_axis_step, v=results[idx][1], 
#             T = results[idx][2], x=results[idx][3],), 
#             DataFrame(results[idx][4], ["term", "level", "state"])
#         )
#     end
#     filename = "/04mult_dynamics_B_δfunction.jld2"
#     save(path*filename, 
#         "df_Result", df_array,
#         "B_list", B_list,
#         "dt", dt,   "t_end", t_end, 
#         "v0", v0,   "δ", δ, "s", s, "ω_p", ω_p,
#         "sampling_rate", sampling_rate,
#         "direction", direction, 
#         "width", width, "shift", shift
#     )
#     return # nothing
# end
# # B のディスパッチ: Function
# function jld_write_mk_single_array(results, dt, t_end, v0, s::Real, δ::Function, B_list::Array{Float64, 1}, 
#     ω_p, width, shift; sampling_rate = 1, direction = "none"::String, path = pwd())

#     println("Export your input and the result by jld file which manipulated δ function")

#     t_axis_step = collect(0:sampling_rate *dt:t_end )
#     df_array= Array{DataFrame}(undef, width)
#     @threads for idx = eachindex(results) # output
#         df_array[idx] = hcat(DataFrame(time=t_axis_step, v=results[idx][1], 
#             T = results[idx][2], x=results[idx][3],), 
#             DataFrame(results[idx][4], ["term", "level", "state"])
#         )
#     end
#     filename = "/04mult_dynamics_B_δfunction.jld2"
#     save(path*filename, 
#         "df_Result", df_array,
#         "B_list", B_list,
#         "dt", dt,   "t_end", t_end, 
#         "v0", v0,   "δ", δ, "s", s, "ω_p", ω_p,
#         "sampling_rate", sampling_rate,
#         "direction", direction, 
#         "width", width, "shift", shift
#     )
#     return # nothing
# end



# # direction のディスパッチ
# function jld_write_mk_single_array(results, dt, t_end, s::Real, δ::Real, 
#     B::Real, ω_p, width, shift; 
#     sampling_rate = 1, direction_list = ["none", "z", "-z"]::Array{String, 1}, path = pwd())

#     println("Export your input and the result by jld file")

#     t_axis_step = collect(0:sampling_rate *dt:t_end )
#     df_array= Array{DataFrame}(undef, width)
#     @threads for idx = eachindex(results) # output
#         df_array[idx] = hcat(DataFrame(time=t_axis_step, v=results[idx][1], 
#             T = results[idx][2], x=results[idx][3],), 
#             DataFrame(results[idx][4], ["term", "level", "state"])
#         )
#     end
#     bgn_str   =    direction_list[begin + shift]
#     terminal_str = direction_list[width + shift]
#     filename = "/05mult_dynamics_direct_"*bgn_str*"-"*terminal_str*".jld2"
#     save(path*filename, 
#         "df_Result", df_array,
#         "direction_list", direction_list, 
#         "dt", dt,   "t_end", t_end, 
#         "s", s, "δ", δ, "B", B, "ω_p", ω_p,
#         "sampling_rate", sampling_rate,
#         "width", width, "shift", shift
#     )
#     return # nothing
# end
# # direction のディスパッチ: Function
# function jld_write_mk_single_array(results, dt, t_end, s::Real, δ::Function, 
#         B::Real, ω_p, width, shift; sampling_rate = 1, 
#         direction_list = ["none", "z", "-z"]::Array{String, 1}, path = pwd())

#     println("Export your input and the result by jld file which manipulated δ function")

#     t_axis_step = collect(0:sampling_rate *dt:t_end )
#     df_array= Array{DataFrame}(undef, width)
#     @threads for idx = eachindex(results) # output
#         df_array[idx] = hcat(DataFrame(time=t_axis_step, v=results[idx][1], 
#             T = results[idx][2], x=results[idx][3],), 
#             DataFrame(results[idx][4], ["term", "level", "state"])
#         )
#     end
#     filename = "/05mult_dynamics_direct_δfunction.jld2"
#     save(path*filename, 
#         "df_Result", df_array,
#         "direction_list", direction_list, 
#         "dt", dt,   "t_end", t_end, 
#         "s", s, "δ", δ, "B", B, "ω_p", ω_p,
#         "sampling_rate", sampling_rate,
#         "width", width, "shift", shift
#     )
#     return # nothing
# end



# # sのディスパッチ
# function jld_write_mk_data_dep_detune(vTx_result, s_list::Array{Float64, 1}, 
#         δ_list::Array{Real, 1}, B::Real, ω_p; start_save_v = t_end, direct="none", path = pwd())
#     println("Export your input and the result of s_list mk_data_dep_detune by jld file")
#     bgn_str   =    string(s_list[begin])
#     terminal_str = string(s_list[end])
#     filename = "/06final_val_depδ_s_"*bgn_str*"-"*terminal_str*".jld2"

#     df_array= Array{DataFrame}(undef, length(s_list))
#     for i = eachindex(s_list)
#         df_array[i] = hcat(DataFrame(δ=δ_list), DataFrame(vTx_result[i], ["EXP_v", "EXP_T", "x_final"]))
#     end
#     save(path*filename, 
#         "df_results", df_array,
#         "δ_list", δ_list, "s_list", s_list, 
#           "B", B, "ω_p", ω_p, "direct", direct, 
#         "record time", start_save_v, 
#     )
#     return # nothing
# end
# # sのディスパッチ: function
# function jld_write_mk_data_dep_detune(vTx_result, s_list::Array{Float64, 1}, 
#         δ_list::Array{Function, 1}, B::Real, ω_p; start_save_v = t_end, direct="none", path = pwd())
    
#         println("Export your input and the result of s_list and δFunctions mk_data_dep_detune by jld file")
#     bgn_str   =    string(s_list[begin])
#     terminal_str = string(s_list[end])
#     filename = "/06final_val_dep_funcδ_s_"*bgn_str*"-"*terminal_str*".jld2"

#     df_array= Array{DataFrame}(undef, length(s_list))
#     for i = eachindex(s_list)
#         df_array[i] = hcat(DataFrame(δ= string.(δ_list)), DataFrame(vTx_result[i], ["EXP_v", "EXP_T", "x_final"]))
#     end
#     save(path*filename, 
#         "df_results", df_array,
#         "δ_list", δ_list, "s_list", s_list, 
#           "B", B, "ω_p", ω_p, "direct", direct, 
#         "record time", start_save_v, 
#     )
#     return # nothing
# end


# # Just Now
# # B のディスパッチ
# function jld_write_mk_data_dep_detune(vTx_result, s::Real, δ_list, 
#         B_list::Array{Float64, 1}, ω_p; start_save_v = t_end, direct="none", path = pwd())
    
#     println("Export your input and the result of B_list mk_data_dep_detune by jld file")
#     bgn_str   =    string(B_list[begin])
#     terminal_str = string(B_list[end])
#     filename = "/07final_val_depδ_B"*bgn_str*"-"*terminal_str*".jld2"

#     df_array= Array{DataFrame}(undef, length(B_list))
#     for i = eachindex(B_list)
#         df_array[i] = hcat(DataFrame(δ=δ_list), DataFrame(vTx_result[i], ["EXP_v", "EXP_T", "x_final"]))
#     end
#     save(path*filename, 
#         "df_results", df_array,
#         "δ_list", δ_list, "B_list", B_list, "ω_p", ω_p, 
#         "s", s, "direct", direct, 
#         "record time", start_save_v
#     )
#     return # nothing
# end
# # B のディスパッチ: function
# function jld_write_mk_data_dep_detune(vTx_result, s::Real, δ_list::Array{Function, 1}, 
#         B_list::Array{Float64, 1}, ω_p; start_save_v = t_end, direct="none", path = pwd())
    
#         println("Export your input and the result of B_list and δFunctions mk_data_dep_detune by jld file")
#     bgn_str   =    string(B_list[begin])
#     terminal_str = string(B_list[end])
#     filename = "/07final_val_dep_funcδ_B"*bgn_str*"-"*terminal_str*".jld2"
    
#     df_array= Array{DataFrame}(undef, length(B_list))
#     for i = eachindex(B_list)
#         df_array[i] = hcat(DataFrame(δ= string.(δ_list)), DataFrame(vTx_result[i], ["EXP_v", "EXP_T", "x_final"]))
#     end
#     save(path*filename, 
#         "df_results", df_array,
#         "δ_list", δ_list, "B_list", B_list, "ω_p", ω_p, 
#         "s", s, "direct", direct, 
#         "record time", start_save_v
#     )
#     return # nothing
# end


# # u0 のディスパッチ
# function jld_write_mk_data_dep_detune(vTx_result, u0_list::Array{Array{Float64, 1}, 1}, s::Real, 
#         δ_list, B::Real, ω_p; start_save_v = t_end, direct="none", path = pwd())

#     println("Export your input and the result of u0_list mk_data_dep_detune by jld file")
#     bgn_str   =    string(u0_list[begin][2])
#     terminal_str = string(u0_list[end][2])
#     filename = "/08final_val_depδ_u0"*bgn_str*"-"*terminal_str*".jld2"

#     df_array = Array{DataFrame}(undef, length(u0_list))
#     for i = eachindex(u0_list)
#         df_array[i] = hcat(DataFrame(δ=δ_list), DataFrame(vTx_result[i], ["EXP_v", "EXP_T", "x_final"]))
#     end
#     save(path*filename, 
#         "df_results", df_array,
#         "δ_list", δ_list, "u0_list", u0_list, 
#         "B", B,  "s", s, "ω_p", ω_p, "direct", direct, 
#         "record time", start_save_v, "sampling_rate", sampling_rate
#     )
#     return # nothing
# end
# # u0 のディスパッチ: function
# function jld_write_mk_data_dep_detune(vTx_result, u0_list::Array{Float64, 1}, s::Real, 
#         δ_list::Array{Function, 1}, B::Real, ω_p; start_save_v = t_end, direct="none", path = pwd())
    
#     println("Export your input and the result of u0_list and δFunctions mk_data_dep_detune by jld file")
#     bgn_str   =    string(u0_list[begin][2])
#     terminal_str = string(u0_list[end][2])
#     filename = "/08final_val_depδ_u0"*bgn_str*"-"*terminal_str*".jld2"

#     df_array = Array{DataFrame}(undef, length(u0_list))
#     for i = eachindex(u0_list)
#         df_array[i] = hcat(DataFrame(δ=string.(δ_list)), DataFrame(vTx_result[i], ["EXP_v", "EXP_T", "x_final"]))
#     end
#     save(path*filename, 
#         "df_results", df_array,
#         "δ_list", δ_list, "u0_list", u0_list, 
#         "B", B, "ω_p", ω_p, "s", s, "direct", direct, 
#         "record time", start_save_v, "sampling_rate", sampling_rate
#     )
#     return # nothing
# end


# # 離調依存性のみ
# function jld_write_mk_data_dep_detune(res_array, s::Real, 
#         δ_list, B::Real, ω_p; start_save_v = t_end, direct="none", path = pwd())

#     println("Export your input and the result of only δ_list mk_data_dep_detune by jld file")
#     bgn_str   =    string(round(δ_list[begin]/(2π*1e6); sigdigits = 2))
#     terminal_str = string(round(δ_list[end]/(2π*1e6); sigdigits = 2))
#     filename = "/09final_val_depδ_"*bgn_str*"-"*terminal_str*".jld2"
#     df = hcat(DataFrame(δ=δ_list), DataFrame(res_array, ["EXP_v", "EXP_T", "x_final"]))
#     save(path*filename, 
#         "df_Result", df,
#         "δ_list", δ_list, 
#         "B", B, "ω_p", ω_p, "s", s, "direct", direct, 
#         "record time", start_save_v
#     )
#     return # nothing
# end


# # 
# function jld_write_eval_trap_prob(exp_vTx_dep_delta, std_dev, prob_results,
#         u0_list, s, δ_list, B, ω_p;
#         start_save_v = t_end, n = 50, direct="none", threshold = 5.0e-3, path = pwd())
#     println("Export your input and the result of eval_trap_prob by jld file")
#     filename = "/10TrapEval_grav_"*direct*".jld2"

#     v0_list = last.(u0_list)

#     df_array = Array{DataFrame}(undef, length(δ_list))
#     for i = eachindex(exp_vTx_dep_delta)
#         df_array[i] = hcat(DataFrame(v0=v0_list), DataFrame(exp_vTx_dep_delta[i], ["EXP_v", "EXP_T", "x_final"]),
#             DataFrame(std_dev[i], ["σ_v", "σ_T", "σ_x"]), 
#             DataFrame(Trap_Probability=prob_results[i])
#         )
#     end
#     save(path*filename, 
#         "df_Result", df_array,
#         "δ_list", δ_list, "u0_list", u0_list, 
#         "B", B, "ω_p", ω_p, "s", s, "direct", direct, 
#         "threshold", threshold,
#         "N", n ,"record time", start_save_v
#     )
#     return # nothing
# end

end # module