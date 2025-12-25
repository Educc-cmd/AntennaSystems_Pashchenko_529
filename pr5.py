import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

lam = 3.0
a = 2.3
N = 6
d = 3.2

print(f"N={N}, d={d} см, λ={lam} см")

lam_cr = 2 * a
if lam < lam_cr:
    lam_xb = lam / np.sqrt(1 - (lam / lam_cr) ** 2)
    print(f"λ_xb = {lam_xb:.4f} см")
else:
    lam_xb = lam
    print("Помилка: λ > λ_cr")

k = 2 * np.pi / lam
beta = 2 * np.pi / lam_xb

theta_deg = np.arange(-90, 90.01, 0.05)
theta_rad = np.radians(theta_deg)

def calc_ArrayFactor(theta):
    Psi = k * d * np.sin(theta) - beta * d
    num = np.sin((N * Psi) / 2)
    den = N * np.sin(Psi / 2)

    val = np.ones_like(theta)
    mask = np.abs(den) > 1e-9
    val[mask] = num[mask] / den[mask]
    return np.abs(val)


def calc_Element_H(theta):
    num = np.cos((np.pi / 2) * np.sin(theta))
    den = np.cos(theta)

    val = np.zeros_like(theta)
    mask = np.abs(den) > 1e-9
    val[mask] = num[mask] / den[mask]
    return np.abs(val)

F_array = calc_ArrayFactor(theta_rad)
F_elem = calc_Element_H(theta_rad)

FE_norm = F_array / np.max(F_array)

FH = F_array * F_elem
FH_norm = FH / np.max(FH)


def find_extrema(angles, values, name):
    max_idxs = argrelextrema(values, np.greater)[0]
    min_idxs = argrelextrema(values, np.less)[0]

    print(f"\n{name}: ТАБЛИЦЯ ЕКСТРЕМУМІВ")
    print(f"{'Тип':<10} | {'Кут (град)':<10} | {'Значення (норм)'}")
    print("-" * 40)

    points_max_x = []
    points_max_y = []
    for i in max_idxs:
        if values[i] > 0.001:
            print(f"{'MAX':<10} | {angles[i]:<10.2f} | {values[i]:.4f}")
            points_max_x.append(angles[i])
            points_max_y.append(values[i])

    print("-" * 40)

    points_min_x = []
    points_min_y = []
    for i in min_idxs:
        print(f"{'MIN':<10} | {angles[i]:<10.2f} | {values[i]:.4f}")
        points_min_x.append(angles[i])
        points_min_y.append(values[i])

    return points_max_x, points_max_y, points_min_x, points_min_y


h_max_x, h_max_y, h_min_x, h_min_y = find_extrema(theta_deg, FH_norm, "ПЛОЩИНА H")
e_max_x, e_max_y, e_min_x, e_min_y = find_extrema(theta_deg, FE_norm, "ПЛОЩИНА E")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

ax1.plot(theta_deg, FH_norm, 'b-', label='FH (теорія)')
ax1.plot(h_max_x, h_max_y, 'ro', markersize=5, label='Максимуми')
ax1.plot(h_min_x, h_min_y, 'ko', markersize=5, label='Мінімуми')
ax1.set_title("Площина H (Вертикальна)")
ax1.set_xlabel('Кут θ, градуси')
ax1.set_ylabel('F(θ)')
ax1.grid(True, linestyle='--', alpha=0.6)
ax1.legend()

ax2.plot(theta_deg, FE_norm, 'g-', label='FE (теорія)')
ax2.plot(e_max_x, e_max_y, 'ro', markersize=5, label='Максимуми')
ax2.plot(e_min_x, e_min_y, 'ko', markersize=5, label='Мінімуми')
ax2.set_title("Площина E (Горизонтальна)")
ax2.set_xlabel('Кут θ, градуси')
ax2.set_ylabel('F(θ)')
ax2.grid(True, linestyle='--', alpha=0.6)
ax2.legend()

plt.tight_layout()
plt.show()