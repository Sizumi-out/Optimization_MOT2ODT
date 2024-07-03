module ExportGraff
    using Plots, LaTeXStrings
    using Base.Threads
    using CSV, DataFrames


    export mk_dynamics_graff
    """
    mk_dynamics_graff(df, dt, t_end, Atom,
            v0::Real, s::Real, δ_list::Array{Float64, 1}, B::Real;\n
            sampling = 1, direct="none"::String, target="position")\n
    ある1つの数について、リストを与えるとそれについて、プロットを重ねる。リストを与えることができるのは、以下の変数。\n
        v0, s, δ, B, direct::string\n
    target: 以下のオプション全てを返す"all"が存在\n
        "position" :(default) 座標のプロットを返す。\n
        "temperture": 温度のプロットを返す。\n
    """
    function mk_dynamics_graff(df_array, dt, t_end, Atom,
            v0::Real, s::Real, δ_list::Array{Float64, 1}, B::Real;
            sampling = 1, direct="none"::String, target="position")
        step_dt = dt * sampling # δのディスパッチ
        t_axis = collect(0:step_dt:t_end - step_dt)
        t_ticks_array = (t_end/4:t_end/4:t_end)
        title_txt = L"s = %$s , B =  %$B, Direction: %$direct, v_0 = %$v0"
        color_list = palette(:darkrainbow, length(δ_list))
        if target == "position" 
            println("make x dynamics")
            plt = plot(leg=:top, xlabel="time[ms]", ylabel="X [m]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = title_txt #, size = (800, 500)
            )
            @threads for idx = eachindex(df_array)
                plot!(plt, t_axis, df_array[idx].x, color=color_list[idx],
                    label="δ = "*string(round(δ_list[idx]/1e6; sigdigits = 3))*"MHz"
                )
            end
        elseif target == "temperture"
            println("make T dynamics")
            plt = plot(leg=:topright, yaxis=:log, xlabel="time[ms]", ylabel="T[K]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = title_txt #, size = (800, 500)
            )
            @threads for idx = eachindex(df_array)
                plot!(plt, t_axis, df_array[idx].T, color=color_list[idx],
                    label="δ = "*string(round(δ_list[idx]/1e6; sigdigits = 3))*"MHz"
                )
            end
        elseif target == "all"
            println("make x and T dynamics")
            plt_posi = plot(leg=:top, xlabel="time[ms]", ylabel="X [m]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = title_txt #, size = (800, 500)
            )
            plt_temp = plot(leg=:topright, xlabel="time[ms]", ylabel="T[K]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = title_txt #, size = (800, 500), yaxis=:log, 
            )
            @threads for idx = eachindex(df_array)
                lbel_txt = "δ = "*string(round(δ_list[idx]/1e6; sigdigits = 3))*"MHz"
                plot!(plt_posi, t_axis, df_array[idx].x, label=lbel_txt, color=color_list[idx])
                plot!(plt_temp, t_axis, df_array[idx].T, label=lbel_txt, color=color_list[idx])
            end
            return plt_posi, plt_temp
        else
            error("no option !!!")
        end
        return plt
    end

    function mk_dynamics_graff(df_array, dt, t_end, Atom,
            v0::Real, s_list::Array{Float64, 1}, δ::Real, B::Real;
            sampling = 1, direct="none"::String, target="position")
        step_dt = dt * sampling # sのディスパッチ
        t_axis = collect(0:step_dt:t_end - step_dt)
        t_ticks_array = (t_end/4:t_end/4:t_end)
        title_txt = L"δ = %$(round(δ/(2*π*1e6); sigdigits = 4))MHz, B = %$B, Direction: %$direct"
        color_list = palette(:darkrainbow, length(s_list))
        if target == "position" 
            println("make x dynamics")
            plt = plot(leg=:top, xlabel="time[ms]", ylabel="X [m]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = title_txt#, size = (800, 500)
            )
            @threads for idx = eachindex(df_array)
                plot!(plt, t_axis, df_array[idx].x, color=color_list[idx],
                    label=L"s = %$(s_list[idx])"
                )
            end
        elseif target == "temperture"
            println("make T dynamics")
            plt = plot(leg=:topright, yaxis=:log, xlabel="time[ms]", ylabel="T[K]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = title_txt #, size = (800, 500)
            )
            @threads for idx = eachindex(df_array)
                plot!(plt, t_axis, df_array[idx].T, color=color_list[idx],
                    label=L"s = %$(s_list[idx])"
                )
            end
        elseif target == "all"
            println("make x and T dynamics")
            plt_posi = plot(leg=:top, xlabel="time[ms]", ylabel="X [m]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = title_txt #, size = (800, 500)
            )
            plt_temp = plot(leg=:topright, xlabel="time[ms]", ylabel="T[K]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = title_txt #, size = (800, 500), yaxis=:log, 
            )
            @threads for idx = eachindex(df_array)
                plot!(plt_posi, t_axis, df_array[idx].x, label=L"s = %$(s_list[idx])", color=color_list[idx])
                plot!(plt_temp, t_axis, df_array[idx].T, label=L"s = %$(s_list[idx])", color=color_list[idx])
            end
            return plt_posi, plt_temp
        else
            error("no option !!!")
        end
        return plt
    end

    # NUMにしたい, Bを[Gauss/cm]にする, NUMにしたい???
    function mk_dynamics_graff(df_array, dt, t_end, Atom,
            v0::Real, s::Real,  δ::Real, B_list::Array{Float64, 1};
            sampling = 1, direct="none"::String, target="position")
        step_dt = dt * sampling # Bのディスパッチ
        t_axis = collect(0:step_dt:t_end - step_dt)
        t_ticks_array = (t_end/4:t_end/4:t_end)
        title_txt = L"δ = %$(round(δ/(2*π*1e6); sigdigits = 4)) MHz, s =  %$s, Direction: %$direct"
        color_list = palette(:darkrainbow, length(B_list))
        if target == "position" 
            println("make x dynamics")
            plt = plot(leg=:top, xlabel="time[ms]", ylabel="X [m]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = title_txt#, size = (800, 500)
            )
            @threads for idx = eachindex(df_array)
                plot!(plt, t_axis, df_array[idx].x, color=color_list[idx],
                    label=L"B = %$(B_list[idx])"
                )
            end
        elseif target == "temperture"
            println("make T dynamics")
            plt = plot(leg=:topright, yaxis=:log, xlabel="time[ms]", ylabel="T[K]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = title_txt #, size = (800, 500)
            )
            @threads for idx = eachindex(df_array)
                plot!(plt, t_axis, df_array[idx].T, color=color_list[idx],
                    label=L"B = %$(B_list[idx])"
                )
            end
        elseif target == "all"
            println("make x and T dynamics")
            plt_posi = plot(leg=:top, xlabel="time[ms]", ylabel="X [m]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = title_txt #, size = (800, 500)
            )
            plt_temp = plot(leg=:topright, xlabel="time[ms]", ylabel="T[K]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = title_txt #, size = (800, 500), yaxis=:log, 
            )
            @threads for idx = eachindex(df_array)
                plot!(plt_posi, t_axis, df_array[idx].x, label=L"B = %$(B_list[idx])", color=color_list[idx])
                plot!(plt_temp, t_axis, df_array[idx].T, label=L"B = %$(B_list[idx])", color=color_list[idx])
            end
            return plt_posi, plt_temp
        else
            error("no option !!!")
        end
        return plt
    end


    function mk_dynamics_graff(df_array, dt, t_end, Atom, v0::Real, s::Real,  δ::Real, B::Real;
        sampling = 1, direct_list=["none"]::Array{String, 1}, target="position")
        step_dt = dt * sampling # Directionのディスパッチ
        t_axis = collect(0:step_dt:t_end - step_dt)
        t_ticks_array = (t_end/4:t_end/4:t_end)
        title_txt = L"δ = %$(round(δ/(2*π*1e6); sigdigits = 4)) MHz, s =  %$s, B: %$B"
        color_list = palette(:darkrainbow, length(direct_list))
        if target == "position" 
            println("make x dynamics")
            plt = plot(leg=:top, xlabel="time[ms]", ylabel="X [m]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = title_txt #, size = (800, 500)
            )
            @threads for idx = eachindex(df_array)
                plot!(plt, t_axis, df_array[idx].x, color=color_list[idx],
                    label="Gravity:"*direct_list[idx]
                )
            end
        elseif target == "temperture"
            println("make T dynamics")
            plt = plot(leg=:topright, yaxis=:log, xlabel="time[ms]", ylabel="T[K]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = title_txt #, size = (800, 500)
            )
            @threads for idx = eachindex(df_array)
                plot!(plt, t_axis, df_array[idx].T, color=color_list[idx],
                    label="Gravity:"*direct_list[idx]
                )
            end
        elseif target == "all"
            println("make x and T dynamics")
            plt_posi = plot(leg=:top, xlabel="time[ms]", ylabel="X [m]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = title_txt #, size = (800, 500)
            )
            plt_temp = plot(leg=:topright, xlabel="time[ms]", ylabel="T[K]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = title_txt #, size = (800, 500), yaxis=:log, 
            )
            @threads for idx = eachindex(df_array)
                plot!(plt_posi, t_axis, df_array[idx].x, label="Gravity:"*direct_list[idx], color=color_list[idx])
                plot!(plt_temp, t_axis, df_array[idx].T, label="Gravity:"*direct_list[idx], color=color_list[idx])
            end
            return plt_posi, plt_temp
        else
            error("no option !!!")
        end
        return plt
    end

    #= 
    function mk_t_graff_single_list(dt, t_end, v0_list::Array{Float64, 1}, s::Real, δ::Real, Atom, AtomState, B::Real; step = 1, start_save_v = t_end, direct="none"::String, target="position")
        step_dt = dt * step # v0のディスパッチ
        t_axis = collect(0:step_dt:t_end - step_dt)
        t_ticks_array = (t_end/4:t_end/4:t_end)
        plt_title = "s = %$s, δ = "*string(round(δ/(2*π*1e6); sigdigits = 4))*"MHz"*", B =  %$B Direction: %$direct"
        results = pmap(v0 -> laser_cooling(dt, t_end, v0, s, δ, Atom, AtomState, B, step = step, v_length = round(Int64, start_save_v/dt), direction = direct), v0_list)
        if target == "position"
            println("make x dynamics")
            plt = plot(leg=:top, xlabel="time[ms]", ylabel="X [m]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = plt_title, size = (800, 500)
            )
            for (i, v0) in enumerate(v0_list)
                plot!(plt, t_axis, results[i][2], label="v0 = %$v0")
            end
        elseif  target == "temperture" 
            println("make T dynamics")
            plt = plot(leg=:topright, yaxis=:log, xlabel="time[ms]", ylabel="T[K]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = plt_title, size = (800, 500)
            )
            for (i, v0) in enumerate(v0_list)
                plot!(plt, t_axis, T_from_v.(results[i][1], (Atom, )), label="v0 = %$v0")
            end
        elseif target == "all"
            println("make x and T dynamics")
            plt_posi = plot(leg=:top, xlabel="time[ms]", ylabel="X [m]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = plt_title, size = (800, 500)
            )
            plt_temp = plot(leg=:topright, yaxis=:log, xlabel="time[ms]", ylabel="T[K]",
                xticks=(t_ticks_array, string.(1e3*t_ticks_array)),
                title = plt_title
            )
            for (i, v0) in enumerate(v0_list)
                plot!(plt_posi, t_axis, results[i][2], label="v0 = %$v0")
                plot!(plt_temp, t_axis, T_from_v.(results[i][1], (Atom, )), label="v0 = %$v0")
            end
            return plt_posi, plt_temp, results
        else
            error("no option !!!")
        end
        return plt, results
    end
    =#
end # end module