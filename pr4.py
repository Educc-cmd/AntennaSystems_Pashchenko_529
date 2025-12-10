import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.special import j1

lambd = 0.03
h = 0.12
l = 0.237
d_cp = 0.023
epsilon_r = 2.2

ksi = 1 + lambd / (2 * l)
print(f"Коефіцієнт уповільнення (ksi): {ksi:.4f}")

k = 2 * np.pi / lambd
print(f"Хвильове число k: {k:.2f}")


def lambda1(x):
    if abs(x) < 1e-9:
        return 1.0
    return 2 * j1(x) / x


steps = []
F_beg_list = []
F1e_list = []
F1h_list = []
Fe_single_list = []
Fh_single_list = []
Fc_list = []
Fe_double_list = []
Fh_double_list = []

for theta_deg in np.arange(0, 90.05, 0.05):
    theta_rad = np.radians(theta_deg)
    sin_theta = np.sin(theta_rad)
    cos_theta = np.cos(theta_rad)

    arg1 = (np.pi * l / lambd) * (ksi - cos_theta)
    term1 = (ksi - 1) / np.sin((np.pi * l / lambd) * (ksi - 1))

    if abs(ksi - cos_theta) < 1e-9:
        val_fb = term1 * (np.pi * l / lambd)
    else:
        val_fb = term1 * np.sin(arg1) / (ksi - cos_theta)
    f_b = abs(val_fb)

    arg_d = (np.pi * d_cp / lambd) * sin_theta
    lam_val = lambda1(arg_d)
    f1e = abs(lam_val * cos_theta)
    f1h = abs(lam_val)

    fe_single = f_b * f1e
    fh_single = f_b * f1h

    arg_c = (np.pi * h / lambd) * sin_theta
    f_c = abs(np.cos(arg_c))

    fh_double = fh_single * f_c
    fe_double = fe_single * f_c

    steps.append(theta_deg)
    F_beg_list.append(f_b)
    F1e_list.append(f1e)
    F1h_list.append(f1h)
    Fe_single_list.append(fe_single)
    Fh_single_list.append(fh_single)
    Fc_list.append(f_c)
    Fe_double_list.append(fe_double)
    Fh_double_list.append(fh_double)


def normalize(data):
    max_val = max(data) if max(data) > 0 else 1
    return [x / max_val for x in data]


F_beg_norm = normalize(F_beg_list)
Fe_s_norm = normalize(Fe_single_list)
Fh_s_norm = normalize(Fh_single_list)
Fe_d_norm = normalize(Fe_double_list)
Fh_d_norm = normalize(Fh_double_list)
Fc_norm = Fc_list


def get_beamwidth(values, angles):
    for i in range(len(values)):
        if values[i] < 0.707:
            y1, y2 = values[i - 1], values[i]
            x1, x2 = angles[i - 1], angles[i]
            angle_0707 = x1 + (0.707 - y1) * (x2 - x1) / (y2 - y1)
            return angle_0707 * 2
    return 0


def find_peaks_min_max(values, angles, threshold=0.05):
    max_x, max_y = [], []
    min_x, min_y = [], []
    passed_main_lobe = False

    for i in range(1, len(values) - 1):
        if not passed_main_lobe:
            if values[i] < values[i - 1] and values[i] < values[i + 1]:
                passed_main_lobe = True
                min_x.append(angles[i])
                min_y.append(0.0)  # Нуль
            continue

        if values[i] > values[i - 1] and values[i] > values[i + 1]:
            max_x.append(angles[i])
            max_y.append(values[i])
        elif values[i] < values[i - 1] and values[i] < values[i + 1]:
            min_x.append(angles[i])
            min_y.append(0.0)

    return min_x, min_y, max_x, max_y

bw_h_s = get_beamwidth(Fh_s_norm, steps)
bw_e_s = get_beamwidth(Fe_s_norm, steps)
min_h_x, min_h_y, max_h_x, max_h_y = find_peaks_min_max(Fh_s_norm, steps)
min_e_x, min_e_y, max_e_x, max_e_y = find_peaks_min_max(Fe_s_norm, steps)

bw_h_d = get_beamwidth(Fh_d_norm, steps)
bw_e_d = get_beamwidth(Fe_d_norm, steps)
min_hd_x, min_hd_y, max_hd_x, max_hd_y = find_peaks_min_max(Fh_d_norm, steps)
min_ed_x, min_ed_y, max_ed_x, max_ed_y = find_peaks_min_max(Fe_d_norm, steps)

print("-" * 40)
print("РЕЗУЛЬТАТИ РОЗРАХУНКУ (ОДНОСТРИЖНЕВА)")
print(f"Ширина ДС (H-площина): {bw_h_s:.2f} град")
print(f"Ширина ДС (E-площина): {bw_e_s:.2f} град")
if len(max_h_y) > 0: print(f"Рівень 1-ї бічної (H): {max_h_y[0]:.3f} (кут {max_h_x[0]:.1f})")
if len(max_e_y) > 0: print(f"Рівень 1-ї бічної (E): {max_e_y[0]:.3f} (кут {max_e_x[0]:.1f})")

print("-" * 40)
print("РЕЗУЛЬТАТИ РОЗРАХУНКУ (ДВОСТРИЖНЕВА)")
print(f"Ширина ДС (H-площина): {bw_h_d:.2f} град")
print(f"Ширина ДС (E-площина): {bw_e_d:.2f} град")
if len(max_hd_y) > 0: print(f"Рівень 1-ї бічної (H): {max_hd_y[0]:.3f} (кут {max_hd_x[0]:.1f})")
if len(max_ed_y) > 0: print(f"Рівень 1-ї бічної (E): {max_ed_y[0]:.3f} (кут {max_ed_x[0]:.1f})")

plt.figure(figsize=(10, 6))
plt.plot(steps, Fh_s_norm, 'k-', linewidth=1.5, label=r'$|F_{H}(\Theta)|$ (Single)')
plt.plot(steps, Fe_s_norm, 'k--', linewidth=1.5, label=r'$|F_{E}(\Theta)|$ (Single)')

for mx, my in zip(max_h_x, max_h_y):
    plt.plot(mx, my, 'ko', markersize=4)
    plt.text(mx, my + 0.02, f"{my:.2f}", ha='center', va='bottom', fontsize=8, color='black')

for mx, my in zip(max_e_x, max_e_y):
    plt.plot(mx, my, 'o', color='gray', markersize=4)
    plt.text(mx, my - 0.04, f"{my:.2f}", ha='center', va='top', fontsize=8, color='gray')

all_mins = sorted(list(set(min_h_x + min_e_x)))
for mx in all_mins:
    plt.plot(mx, 0, 'bo', markersize=4, label='_nolegend_')

plt.axhline(y=0.707, color='gray', linestyle=':', linewidth=0.8)
plt.text(2, 0.72, f"Wid(H)={bw_h_s:.1f}$^\circ$", fontsize=9)
plt.text(2, 0.65, f"Wid(E)={bw_e_s:.1f}$^\circ$", fontsize=9)

plt.title(f"Рис. 1. ДС Однострижневої антени (H та E)\nВаріант 15: $\lambda={lambd * 100}$см, $L={l * 100}$см")
plt.xlabel(r'$\Theta^\circ$')
plt.ylabel('Нормована амплітуда')
plt.xlim(0, 90)
plt.ylim(0, 1.1)
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend([r'$|F_{H}(\Theta)|$ (Single)', r'$|F_{E}(\Theta)|$ (Single)'],
           loc="upper right", bbox_to_anchor=(1, 1), fontsize=10)
plt.savefig("DS_Single_Rod.jpg", dpi=600)
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(steps, Fh_s_norm, 'k-', linewidth=1.2, label=r'$|F_{elem}(\Theta)|$ (Single H)')
plt.plot(steps, Fc_norm, 'k-', linewidth=1.2, label=r'$|F_{C}(\Theta)|$ (Array Factor)')
plt.plot(steps, Fh_d_norm, 'k--', linewidth=2, label=r'$|F_{H}(\Theta)|$ (Total Double)')

for mx, my in zip(max_hd_x, max_hd_y):
    plt.plot(mx, my, 'ko', markersize=4)
    plt.text(mx, my + 0.03, f"{my:.2f}", ha='center', fontsize=8)

for mx, my in zip(min_hd_x, min_hd_y):
    plt.plot(mx, my, 'bo', markersize=4, markerfacecolor='white')

plt.axhline(y=0.707, color='gray', linestyle=':', linewidth=0.8)
plt.text(2, 0.72, f"Wid={bw_h_d:.1f}$^\circ$", fontsize=9)

plt.title(f"Рис. 2. ДС Двострижневої антени (Площина H)\nВаріант 15: $h={h * 100}$см")
plt.xlabel(r'$\Theta^\circ$')
plt.ylabel('Нормована амплітуда')
plt.xlim(0, 90)
plt.ylim(0, 1.1)
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend(loc="upper right", bbox_to_anchor=(1, 1), fontsize=10)
plt.savefig("DS_Double_Rod_H.jpg", dpi=600)
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(steps, Fe_s_norm, 'k-', linewidth=1.2, label=r'$|F_{elem}(\Theta)|$ (Single E)')
plt.plot(steps, Fc_norm, 'k-', linewidth=1.2, label=r'$|F_{C}(\Theta)|$ (Array Factor)')
plt.plot(steps, Fe_d_norm, 'k--', linewidth=2, label=r'$|F_{E}(\Theta)|$ (Total Double)')

for mx, my in zip(max_ed_x, max_ed_y):
    plt.plot(mx, my, 'ko', markersize=4)
    plt.text(mx, my + 0.03, f"{my:.2f}", ha='center', fontsize=8)

for mx, my in zip(min_ed_x, min_ed_y):
    plt.plot(mx, my, 'bo', markersize=4, markerfacecolor='white')

plt.axhline(y=0.707, color='gray', linestyle=':', linewidth=0.8)
plt.text(2, 0.72, f"Wid={bw_e_d:.1f}$^\circ$", fontsize=9)

plt.title(f"Рис. 3. ДС Двострижневої антени (Площина E)\nВаріант 15: $h={h * 100}$см")
plt.xlabel(r'$\Theta^\circ$')
plt.ylabel('Нормована амплітуда')
plt.xlim(0, 90)
plt.ylim(0, 1.1)
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend(loc="upper right", bbox_to_anchor=(1, 1), fontsize=10)
plt.savefig("DS_Double_Rod_E.jpg", dpi=600)
plt.show()