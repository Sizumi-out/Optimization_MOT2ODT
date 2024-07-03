"""運動方程式を連立で解く
u[1] = x, u[2] = v
ẋ = v, v̇ = -1/m * dU/dx
"""
simul_diff_eq(t, u, lase1064, Atom, tθs) = [u[2], -tθs.acc - grad_Us(t, u[1], lase1064, Atom, tθs)/Atom.m] # better

get_gravity(tgs::tgFn, cnt) =  gravity(tgs.g_sw_ary[cnt])
get_gravity(tgs::tg, cnt) = tgs.g

"""rk_update!(tvx, tg, lase1064, Atom): Fourth-order Runge-Kutta method"""
function rk_update!(tvx, tθs, tg, Atom, lase1064)
    tθs.acc = get_gravity(tg, tvx.cnt)
    tθs.t = tg.dt*(tvx.cnt - 1)
    du1 = simul_diff_eq(tθs.t,              tvx.u                 , lase1064, Atom, tθs)
    du2 = simul_diff_eq(tθs.t + tg.dt_half, tvx.u + du1*tg.dt_half, lase1064, Atom, tθs)
    tθs.t += tg.dt_half # to avoid integration
    du3 = simul_diff_eq(tθs.t,              tvx.u + du2*tg.dt_half, lase1064, Atom, tθs)
    du4 = simul_diff_eq(tθs.t + tg.dt_half, tvx.u + du3*tg.dt     , lase1064, Atom, tθs)
    tvx.u += (du1 + 2du2 + 2du3 + du4)*tg.dt/6
    return # nothing
end


function uni_lin_acc!(tvx, tg)
    gr = get_gravity(tg, tvx.cnt)
    @views tvx.u += tg.dt * SVector(tvx.u[2] + gr*tg.dt_half, gr) # [x, v] += [v + g*dt²/2, g*dt]
    return #
end