module Constants

export ħ, k_b, μ_B, gra_acc, c
export ω_from_λ, λ_from_ω, k_from_λ

const ħ = 6.62607004e-34/(2*π)  # プランク定数
const k_b = 1.38064852e-23      # ボルツマン定数
const μ_B = 1.4e6               # [MHz/Gauss] Bohr
const gra_acc = 9.80665         # [m/s^2] 重力加速度
const c = 299_792_458           # 光速

ω_from_λ(λ) = 2*pi*c/λ

λ_from_ω(ω) = 2*pi*c/ω

k_from_λ(λ) = 2*pi/λ

end # end module