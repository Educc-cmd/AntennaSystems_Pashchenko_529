import numpy
import matplotlib.pyplot as plt
import math
from scipy.special import jv

lambd = 0.03
D = 0.8
f = 0.3

print(f"λ = {lambd} (m)")
print(f"D = {D} (m)")
print(f"f = {f} (m)")

k = (2 * numpy.pi) / lambd
R0 = D / 2
p = 2 * f

psi0_rad = 2 * numpy.arctan(R0 / p)
psi0_deg = numpy.degrees(psi0_rad)
cos_psi0 = numpy.cos(psi0_rad)

print(f"R0 = {R0} (m)")
print(f"p = {p} (m)")
print(f"ψ0 = {psi0_deg:.2f}°")

Th = ((1 + cos_psi0) / 2) * numpy.sin((numpy.pi / 2) * cos_psi0)
Te = Th * cos_psi0

print(f"П'єдестал Th (площина H) = {Th:.4f}")
print(f"П'єдестал Te (площина E) = {Te:.4f}")

def Lambda1(u):
    if u == 0: return 1.0
    return (2 * jv(1, u)) / u


def Lambda2(u):
    if u == 0: return 1.0
    return (8 * jv(2, u)) / (u ** 2)

FH_vals = []
FE_vals = []
steps = []

SGP_H = 0
SGP_E = 0
val_SGP_H = 0
val_SGP_E = 0

for teta in numpy.arange(0, numpy.radians(90), 0.0001):

    u = k * R0 * numpy.sin(teta)

    elem_factor = (numpy.cos(teta / 2)) ** 2

    bracket_H = Th * Lambda1(u) + ((1 - Th) / 2) * Lambda2(u)
    val_H = (2 / (1 + Th)) * bracket_H * elem_factor

    bracket_E = Te * Lambda1(u) + ((1 - Te) / 2) * Lambda2(u)
    val_E = (2 / (1 + Te)) * bracket_E * elem_factor

    FH_vals.append(val_H)
    FE_vals.append(val_E)
    deg_teta = numpy.degrees(teta)
    steps.append(deg_teta)

    target_level = 0.707

    if len(FH_vals) > 1:
        if FH_vals[-2] >= target_level and FH_vals[-1] < target_level:
            SGP_H = 2 * deg_teta
            val_SGP_H = val_H

    if len(FE_vals) > 1:
        if FE_vals[-2] >= target_level and FE_vals[-1] < target_level:
            SGP_E = 2 * deg_teta
            val_SGP_E = val_E

print("Ширина головної пелюстки в площині H = " + str(round(SGP_H, 2)) + '\u00b0')
print("Ширина головної пелюстки в площині E = " + str(round(SGP_E, 2)) + '\u00b0')

max_x_FH = []
max_y_FH = []
max_x_FE = []
max_y_FE = []

min_x_FH = []
min_y_FH = []
min_x_FE = []
min_y_FE = []

for i in range(1, len(FH_vals) - 1):
    if FH_vals[i] > FH_vals[i - 1] and FH_vals[i] > FH_vals[i + 1]:
        if steps[i] > 0.5:
            max_x_FH.append(steps[i])
            max_y_FH.append(FH_vals[i])
    if abs(FH_vals[i]) < abs(FH_vals[i - 1]) and abs(FH_vals[i]) < abs(FH_vals[i + 1]):
        min_x_FH.append(steps[i])
        min_y_FH.append(FH_vals[i])

for i in range(1, len(FE_vals) - 1):
    if FE_vals[i] > FE_vals[i - 1] and FE_vals[i] > FE_vals[i + 1]:
        if steps[i] > 0.5:
            max_x_FE.append(steps[i])
            max_y_FE.append(FE_vals[i])
    if abs(FE_vals[i]) < abs(FE_vals[i - 1]) and abs(FE_vals[i]) < abs(FE_vals[i + 1]):
        min_x_FE.append(steps[i])
        min_y_FE.append(FE_vals[i])

len_value = [len(max_x_FH), len(max_x_FE), len(min_x_FH), len(min_x_FE)]
print("\nТабл. 1 - Аналіз ДС Дзеркальної антени в площині Н та Е")
print("-" * 45)
print("| № |θminH|θminE|θmaxH|FH(θ)|θmaxE|FE(θ)|")
for i in range(0, max(len_value)):
    v1 = f" {i + 1:.0f}"
    if len(v1) < len(str(max(len_value))) + 2: v1 += " "

    v2 = f"{min_x_FH[i]:.2f}" if i < len_value[2] else '  -  '
    if len(v2) < 5: v2 += "0"
    v3 = f"{min_x_FE[i]:.2f}" if i < len_value[3] else '  -  '
    if len(v3) < 5: v3 += "0"

    v4 = f"{max_x_FH[i]:.2f}" if i < len_value[0] else '  -  '
    if len(v4) < 5: v4 += "0"
    h4 = f"{max_y_FH[i]:.3f}" if i < len_value[0] else '  -  '
    if len(h4) < 5: h4 += "0"

    v5 = f"{max_x_FE[i]:.2f}" if i < len_value[1] else '  -  '
    if len(v5) < 5: v5 += "0"
    h5 = f"{max_y_FE[i]:.3f}" if i < len_value[1] else '  -  '
    if len(h5) < 5: h5 += "0"

    print(f"|{v1}|{v2}|{v3}|{v4}|{h4}|{v5}|{h5}|")
print("-" * 45)

fig, ax = plt.subplots(figsize=(20 / 2.54, 12 / 2.54))

ax.plot(steps, FH_vals, linewidth=1.2, color='blue', label="$ F_{H}(θ) $ - Площина H")
ax.plot(steps, FE_vals, linewidth=1.2, color='red', linestyle='--', label="$ F_{E}(θ) $ - Площина E")

ax.plot(SGP_H / 2, val_SGP_H, 'bo', markersize=4)
ax.plot(SGP_E / 2, val_SGP_E, 'ro', markersize=4)

plt.annotate(f'H: {SGP_H:.2f}\u00b0', xy=(SGP_H / 2, val_SGP_H), xytext=((SGP_H / 2) + 0.5, val_SGP_H + 0.05),
             arrowprops=dict(arrowstyle='->', color='blue'), fontsize=8, color='blue')
plt.annotate(f'E: {SGP_E:.2f}\u00b0', xy=(SGP_E / 2, val_SGP_E), xytext=((SGP_E / 2) + 0.5, val_SGP_E - 0.1),
             arrowprops=dict(arrowstyle='->', color='red'), fontsize=8, color='red')

ax.hlines(y=0.707, xmin=0, xmax=max(steps), colors='green', linestyles='dotted', linewidth=0.8, label="Рівень 0.707")

ax.plot(max_x_FH, max_y_FH, "o", markersize=5, color="blue", label="$ \\theta_{max} $ FH")
ax.plot(max_x_FE, max_y_FE, "o", markersize=5, color="red", label="$ \\theta_{max} $ FE")
ax.plot(min_x_FE, min_y_FE, "o", markersize=5, color="green", label="$ \\theta_{min} $ FE")
ax.plot(min_x_FH, min_y_FH, "o", markersize=5, color="grey", label="$ \\theta_{min} $ FH")

ax.set_xlabel('Кут θ (градуси)', fontsize=10)
ax.set_ylabel('Нормоване значення поля |F(θ)|', fontsize=10)
ax.set_title(f'ДС Дзеркальної антени (Варіант 15)\nD={D}м, f={f}м, λ={lambd}м', fontsize=11)

plt.xticks(numpy.arange(0, 91, 4), fontsize=8)
plt.yticks(numpy.arange(-0.3, 1.1, 0.1), fontsize=8)
plt.ylim(-0.22, 1.02)
plt.xlim(0, 90)

plt.legend(loc="upper right", fontsize=8)
plt.grid(which='both', linestyle='--', linewidth=0.5, color='gray', alpha=0.7)

fig.savefig("ДС_Дзеркальна_антена_Вар15.jpg", dpi=300)
plt.show()