progress(cnt, axis_len) = round(Int, 1000*cnt/axis_len)


processTrap!(ctrler::fCtrler, tvx, tθs, tg, Atom, state, TrapStrs) = applyTrap!(tvx, tθs, tg, Atom, state, TrapStrs[ctrler.swMat[:, tvx.cnt]])
processTrap!(ctrler::Ctrler, tvx, tθs, tg, Atom, state, TrapStrs) = applyTrap!(tvx, tθs, tg, Atom, state, TrapStrs)


function applyTrap!(tvx, tθs, tg, Atom, AtomState, applyTrapStrs)
    cnt = Int8(1)
    for TrapStr in applyTrapStrs
        # println("@applyTrap! t=",tvx.cnt, "  TrapStr: ", TrapStr, "\t type:", typeof(TrapStr))
        if Trapping!(TrapStr, tvx, tθs, tg, Atom, AtomState, applyTrapStrs, cnt)
            return false
        end
        cnt += 1
    end
    return true
end

function Trapping!(::Dissipative, tvx, tθs, tg, Atom, AtomState, TrapStrs, cnt) 
    MOT!(tvx, tg.dt, TrapStrs[cnt], Atom, AtomState)
    return false
end
function Trapping!(::FORT, tvx, tθs, tg, Atom, AtomState, TrapStrs, cnt) 
    @views rk_update!(tvx, tθs, tg, Atom, TrapStrs[cnt:end]) # todo: @views is effective?
    return true
end


function trap!(tvxs, tg, ctrler, recorder, Atom, trapStructs)
    
    # initial
    tθs = mktθ(trapStructs, tg)

    state = AtomState(1, 1//2, -1//2)
    judge_record(tvxs, ctrler) && record!(recorder, tvxs, state) # record initial states
    tvxs.cnt += 1

    ctrler.prog_sw && (prog = ProgressUnknown("Progress(‰):" ))
    while tvxs.cnt <= ctrler.t_ax_len
        ctrler.prog_sw && ProgressMeter.update!(prog, progress(tvxs.cnt, ctrler.t_ax_len))
        
        # if only CMOT: x update! by gravity & updated v by MOT
        if processTrap!(ctrler, tvxs, tθs, tg, Atom, state, trapStructs)
            uni_lin_acc!(tvxs, tg) # or apply `rk_update!` with lase1064 = 0
        end

        judge_record(tvxs, ctrler) && record!(recorder, tvxs, state)
        tvxs.cnt += 1
    end
    ctrler.prog_sw && ProgressMeter.finish!(prog)
    # println("integral: ", tθs.θ_ω)
    return recorder
end