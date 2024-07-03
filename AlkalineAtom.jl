using LaTeXStrings

export γ_trns_streng_abs, γ_trns_streng_down, state_list
export q_state_cnt_from_array

## Atom select
function select_atom(str::String)
    if str == "6Li"
        return Atom(Li_mass, Li_λ, Li_γ, γ_trns_streng_abs, γ_trns_streng_down)
    elseif  str == "23Na"
        return Atom(Na_mass, Na_λ, Na_γ, γ_trns_streng_abs, γ_trns_streng_down)
    else
        throw(DomainError("There is no atom"))
    end
end

select_atom(x::Atom) = x


## Atoms
# 6Li
const Li_λ = 671 * 1e-9 # 670 [nm]
const Li_mass = 6e-3 / 6.02 *1e-23 # 質量数6[g/mol]
const Li_γ = 2*π* 6e6 # [/rad]


# 23Na
const Na_λ = 589 * 1e-9 # 589 [nm]
const Na_mass = 23e-3 / 6.02 *1e-23 # 質量数23[g/mol]
const Na_γ = 2*π/(16e-9) # [/rad]


## states: clebsh gordan coeeficients

const γ_trns_streng_abs = [
    ([160, 50,  0], [0, 150,   0]),# mF_1per4inside, mF_1per4outside
    ([  5, 64, 81], [0,  48, 162]),# mF_3per4inside, mF_3per4outside
    ([ 15, 48, 27], [0,   0, 270]) # mF_9per4inside, mF_9per4outside
]

const γ_trns_streng_down = [
    [ 80, 160,  10,  5,  15], # F = 1/2, m = ±1/2から
    [100,  50,   8, 64,  48], # F = 3/2, m = ±1/2から
    [150,   0,  48,  0,  72], # F = 3/2, m = ±3/2から
    [  0,   0, 162, 81,  27], # F = 5/2, m = ±1/2から
    [  0,   0, 162,  0, 108], # F = 5/2, m = ±3/2から
    [  0,   0, 270,  0,   0], # F = 5/2, m = ±5/2から
]

function q_state_cnt_from_array(qs_array)
    qs_cnt_array = zeros(Int, 18)
    for i in 1:size(qs_array)[1]
        if qs_array[i, 1] == 1
            (qs_array[i, 2] == 1//2 && qs_array[i, 3] == 1//2) &&(qs_cnt_array[1] += 1)
            (qs_array[i, 2] == 1//2 && qs_array[i, 3] == -1//2) &&(qs_cnt_array[2] += 1)
            (qs_array[i, 2] == 3//2 && qs_array[i, 3] == 1//2) &&(qs_cnt_array[3] += 1)
            (qs_array[i, 2] == 3//2 && qs_array[i, 3] == -1//2) &&(qs_cnt_array[4] += 1)
            (qs_array[i, 2] == 3//2 && qs_array[i, 3] == 3//2) &&(qs_cnt_array[5] += 1)
            (qs_array[i, 2] == 3//2 && qs_array[i, 3] == -3//2) &&(qs_cnt_array[6] += 1)
        elseif qs_array[i, 1] == 2
            (qs_array[i, 2] == 1//2 && qs_array[i, 3] == 1//2) && (qs_cnt_array[7] += 1)
            (qs_array[i, 2] == 1//2 && qs_array[i, 3] == -1//2) && (qs_cnt_array[8] += 1)
            (qs_array[i, 2] == 3//2 && qs_array[i, 3] == 1//2) && (qs_cnt_array[9] += 1)
            (qs_array[i, 2] == 3//2 && qs_array[i, 3] == -1//2) && (qs_cnt_array[10] += 1)
            (qs_array[i, 2] == 3//2 && qs_array[i, 3] == 3//2) && (qs_cnt_array[11] += 1)
            (qs_array[i, 2] == 3//2 && qs_array[i, 3] == -3//2) && (qs_cnt_array[12] += 1)
            (qs_array[i, 2] == 5//2 && qs_array[i, 3] == 1//2) && (qs_cnt_array[13] += 1)
            (qs_array[i, 2] == 5//2 && qs_array[i, 3] == -1//2) && (qs_cnt_array[14] += 1)
            (qs_array[i, 2] == 5//2 && qs_array[i, 3] == 3//2) && (qs_cnt_array[15] += 1)
            (qs_array[i, 2] == 5//2 && qs_array[i, 3] == -3//2) && (qs_cnt_array[16] += 1)
            (qs_array[i, 2] == 5//2 && qs_array[i, 3] == 5//2) && (qs_cnt_array[17] += 1)
            (qs_array[i, 2] == 5//2 && qs_array[i, 3] == -5//2) && (qs_cnt_array[18] += 1)
        else
            error("This result has nonexistent states!!")
        end
    end
    return qs_cnt_array
end

const state_list =[
    L"^2 S_{1/2}, F=1/2, m=1/2", 
    L"^2 S_{1/2}, F=1/2, m=-1/2",
    L"^2 S_{1/2}, F=3/2, m=1/2", 
    L"^2 S_{1/2}, F=3/2, m=-1/2",
    L"^2 S_{1/2}, F=3/2, m=3/2", 
    L"^2 S_{1/2}, F=3/2, m=-3/2",
    L"^2 P_{3/2}, F=1/2, m=1/2", 
    L"^2 P_{3/2}, F=1/2, m=-1/2", 
    L"^2 P_{3/2}, F=3/2, m=1/2", 
    L"^2 P_{3/2}, F=3/2, m=-1/2", 
    L"^2 P_{3/2}, F=3/2, m=3/2", 
    L"^2 P_{3/2}, F=3/2, m=-3/2",
    L"^2 P_{3/2}, F=5/2, m=1/2", 
    L"^2 P_{3/2}, F=5/2, m=-1/2", 
    L"^2 P_{3/2}, F=5/2, m=3/2", 
    L"^2 P_{3/2}, F=5/2, m=-3/2", 
    L"^2 P_{3/2}, F=5/2, m=5/2", 
    L"^2 P_{3/2}, F=5/2, m=-5/2",
]
