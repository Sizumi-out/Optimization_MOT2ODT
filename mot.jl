## MOT functions
# MOT
const g_F = [2.0, 4/3]

# rate calculation
δ_mot(tvx, t, mot::MOT, AtomState, pol) = -(mot.δ + 2*pi * μ_B *((AtomState.m_spin + pol) * g_F[2] - AtomState.m_spin * g_F[1] )* mot.B * tvx.u[1])
γ_abs(tvx, t, mot::MOT, Atom, AtomState, pol) = Atom.γ/(2*2*pi) *mot.s/(1 + mot.s + (2*(δ_mot(tvx, t, mot, AtomState, pol) - pol*Atom.k * tvx.u[2])/Atom.γ)^2 ) # 吸収レートの関数 L: σ+, R: σ- :(2*2*pi): γ[/rad]→[Hz]に補正するために2piで割る


function δ_mot(tvx, t, mot::fMOT, AtomState, pol) 
    δ = get_val_t(mot, :δ, t) # δ = mot.δ(t)
    B = get_val_t(mot, :B, t) # B = mot.B(t)
    return -(δ + 2*pi * μ_B *((AtomState.m_spin + pol) * g_F[2] - AtomState.m_spin * g_F[1] )* B * tvx.u[1])
end
function γ_abs(tvx, dt, mot::fMOT, Atom, AtomState, pol)
    t_shift = dt*(tvx.cnt - 1) - mot.sect[1]
    s = get_val_t(mot, :s, t_shift) # s = mot.s(t_shift)
    return Atom.γ/(2*2*pi) *s/(1 + s + (2*(δ_mot(tvx, t_shift, mot, AtomState, pol) - pol*Atom.k * tvx.u[2])/Atom.γ)^2 )
end


# MOT processing
function absorb!(tvx, dt, mot, Atom, state) # todo: struct or current? (σplus, γ_abs_idx, γ_abs_tpl, norm_val)
    absorb_cond = rand()
    γ_abs_idx = abs_idx(state)
    
    γ_abs_tpl = abs_tpl_idx(state, σplus) # σ+とスピンの情報よりTuple内(1or2)を指定
    norm_val = sum_strenth(Atom, γ_abs_idx, γ_abs_tpl)
    if 0 <= absorb_cond < sum_strenth(Atom, γ_abs_idx, γ_abs_tpl)/sumγ* γ_abs(tvx, dt, mot, Atom, state, σplus) * dt
        abs_atom_state!(state, Atom, σplus, γ_abs_idx, γ_abs_tpl, norm_val) # σ+吸収に伴うstateの更新
        tvx.u[2] += Atom.v_k # 運動量キック
        return # Int8(2) # count for absorption of the σ+
    end

    γ_abs_tpl = abs_tpl_idx(state, σmns) # σ-とスピンの情報よりTuple内(1or2)を指定
    norm_val = sum_strenth(Atom, γ_abs_idx, γ_abs_tpl)
    if 1 - sum_strenth(Atom, γ_abs_idx, γ_abs_tpl)/sumγ*γ_abs(tvx, dt, mot, Atom, state, σmns) * dt < absorb_cond <= 1
        abs_atom_state!(state, Atom, σmns, γ_abs_idx, γ_abs_tpl, norm_val) # σ-吸収に伴うstateの更新
        tvx.u[2] -= Atom.v_k # 運動量キック
    end
    return # Int8(1) # count for absorption of the σ-
end

function emit!(tvx, dt, Atom, state)
    if rand() < Atom.γ*dt/(2*π)
        down_atom_state!(state, Atom)
        tvx.u[2] += Atom.v_k * cos(2*π*rand()) # 遷移(2→1)分の運動量キック
    end
    return 
end

function MOT!(tvx, dt, mot, Atom, state)
    if state.term == 1//1
        absorb!(tvx, dt, mot, Atom, state)
    else
        emit!(tvx, dt, Atom, state)
    end
    return #
end