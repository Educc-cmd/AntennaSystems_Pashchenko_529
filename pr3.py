import numpy
import matplotlib.pyplot as plt
import math
from scipy.signal import argrelextrema

print("\nРОЗРАХУНОК ДЛЯ ПЛОЩИНИ E")

F1_list = [1]
FC_list = [1]
FE_list = [1]
steps = [0]

SGP_FC = 0
SGP_FE = 0
fS_FC = 0
fS_FE = 0

max_x_FC = []
max_y_FC = []
max_x_FE = []
max_y_FE = []

min_x_FC = []
min_y_FC = []
min_x_FE = []
min_y_FE = []

lambd = 0.03
ap = 0.12
bp = 0.14
print("λ = " + str(lambd) + "(m)")
print("bp = " + str(bp) + "(m)")
k = (2 * numpy.pi) / lambd

for teta in numpy.arange(0.0001, numpy.pi / 2, 0.0001):
    mn1 = abs((1 + numpy.cos(teta)) / 2)

    arg = (numpy.pi * bp / lambd) * numpy.sin(teta)
    mn2 = abs(numpy.sin(arg) / arg)

    mn3 = mn1 * mn2

    F1_list += [mn1]
    FC_list += [mn2]
    FE_list += [mn3]

    if 0.7065 < mn2 < 0.7075:
        SGP_FC = 2 * math.degrees(teta)
        fS_FC = mn2

    if 0.7065 < mn3 < 0.7075:
        SGP_FE = 2 * math.degrees(teta)
        fS_FE = mn3

    steps += [math.degrees(teta)]

for i in range(1, len(FC_list) - 1):
    if FC_list[i] > FC_list[i - 1] and FC_list[i] > FC_list[i + 1]:
        max_x_FC.append(steps[i])
        max_y_FC.append(FC_list[i])
    if FE_list[i] > FE_list[i - 1] and FE_list[i] > FE_list[i + 1]:
        max_x_FE.append(steps[i])
        max_y_FE.append(FE_list[i])

for i in range(1, len(FC_list) - 1):
    if FC_list[i] < FC_list[i - 1] and FC_list[i] < FC_list[i + 1]:
        min_x_FC.append(steps[i])
        min_y_FC.append(FC_list[i])
    if FE_list[i] < FE_list[i - 1] and FE_list[i] < FE_list[i + 1]:
        min_x_FE.append(steps[i])
        min_y_FE.append(FE_list[i])

print("Ширина головної пелюстки в площині E = " + str(round(SGP_FE, 2)) + '\u00b0')

len_value = [len(max_x_FC), len(max_x_FE), len(min_x_FC), len(min_x_FE)]
print("Табл. 1 - Аналіз ДС Рупора в площині Е")
print("-----------------------------------------")
print("| № |θminC|θminE|θmaxC|FC(θ)|θmaxE|FE(θ)|")
for i in range(0, max(len_value)):
    v1 = f" {i + 1:.0f}"
    if len(v1) < len(str(max(len_value))) + 2: v1 += " "
    v2 = f"{min_x_FC[i]:.2f}" if i < len(min_x_FC) else '  -  '
    if len(v2) < 5: v2 += "0"
    v3 = f"{min_x_FE[i]:.2f}" if i < len(min_x_FE) else '  -  '
    if len(v3) < 5: v3 += "0"
    v4 = f"{max_x_FC[i]:.2f}" if i < len(max_x_FC) else '  -  '
    if len(v4) < 5: v4 += "0"
    h4 = f"{max_y_FC[i]:.3f}" if i < len(max_x_FC) else '  -  '
    if len(h4) < 5: h4 += "0"
    v5 = f"{max_x_FE[i]:.2f}" if i < len(max_x_FE) else '  -  '
    if len(v5) < 5: v5 += "0"
    h5 = f"{max_y_FE[i]:.3f}" if i < len(max_x_FE) else '  -  '
    if len(h5) < 5: h5 += "0"
    print(f"|{v1}|{v2}|{v3}|{v4}|{h4}|{v5}|{h5}|")
print("-----------------------------------------")

fig1, ax1 = plt.subplots(figsize=(20 / 2.54, 12 / 2.54))
ax1.plot(steps, F1_list, linewidth=0.7, label=r"$ F_{1e}(\theta) $ (Element)")
ax1.plot(steps, FC_list, linewidth=0.7, label=r"$ F_{CE}(\theta) $ (Aperture)")
ax1.plot(steps, FE_list, linewidth=0.7, label=r"$ F_{E}(\theta) $ (Total)")
ax1.plot(SGP_FE / 2, fS_FE, 'go', markersize=4, label="ШГП в площині E")

ax1.annotate(f'({fS_FE:.3f}, {SGP_FE / 2:.2f}\u00b0)',
             xy=(SGP_FE / 2, fS_FE),
             xytext=((SGP_FE / 2) + 3, fS_FE + 0.05),
             arrowprops=dict(arrowstyle='->', color='black'), fontsize=8)

ax1.plot([0, SGP_FE / 2], [fS_FE, fS_FE], 'g--', linewidth=0.5)
ax1.plot([SGP_FE / 2, SGP_FE / 2], [0, fS_FE], 'g--', linewidth=0.5)

ax1.plot(max_x_FC, max_y_FC, "o", markersize=4, color="black", label=r"$ \theta_{max} $ FCE")
ax1.plot(max_x_FE, max_y_FE, "o", markersize=4, color="grey", label=r"$ \theta_{max} $ FE")
ax1.plot(min_x_FE, min_y_FE, "o", markersize=4, color="blue", label=r"$ \theta_{min} $ FE")

ax1.set_xlabel('θ' + '\u00b0', fontsize=10)
ax1.set_ylabel(r'|Fe($\theta$)|', fontsize=10)
plt.xticks(numpy.arange(0, 100, 2), fontsize=7)
plt.yticks(numpy.arange(0, 1.2, 0.1), fontsize=7)
plt.ylim(-0.01, 1.01)
plt.xlim(0, 90.5)
plt.title("Діаграма спрямованості в площині E")
plt.legend(loc="upper right", fontsize=7)
plt.grid(which='both', linestyle='--', linewidth=0.2, color='gray')
fig1.savefig("DS_Plane_E.jpg", dpi=600)

print("\nРОЗРАХУНОК ДЛЯ ПЛОЩИНИ H")

F1_list = [1]
FC_list = [1]
FH_list = [1]
steps = [0]

SGP_FC = 0
SGP_FH = 0
fS_FC = 0
fS_FH = 0

max_x_FC = []
max_y_FC = []
max_x_FH = []
max_y_FH = []

min_x_FC = []
min_y_FC = []
min_x_FH = []
min_y_FH = []

print("ap = " + str(ap) + "(m)")

for teta in numpy.arange(0.0001, numpy.pi / 2, 0.0001):
    mn1 = abs((1 + numpy.cos(teta)) / 2)

    arg = (numpy.pi * ap / lambd) * numpy.sin(teta)
    denom = 1 - ((2 * ap * numpy.sin(teta)) / lambd) ** 2

    if abs(denom) < 0.0001:
        mn2 = numpy.pi / 4
    else:
        mn2 = abs(numpy.cos(arg) / denom)

    mn3 = mn1 * mn2

    F1_list += [mn1]
    FC_list += [mn2]
    FH_list += [mn3]

    if 0.7065 < mn2 < 0.7075:
        SGP_FC = 2 * math.degrees(teta)
        fS_FC = mn2

    if 0.7065 < mn3 < 0.7075:
        SGP_FH = 2 * math.degrees(teta)
        fS_FH = mn3

    steps += [math.degrees(teta)]

for i in range(1, len(FC_list) - 1):
    if FC_list[i] > FC_list[i - 1] and FC_list[i] > FC_list[i + 1] and steps[i] > 5:  # ігноруємо шум біля 0
        max_x_FC.append(steps[i])
        max_y_FC.append(FC_list[i])
    if FH_list[i] > FH_list[i - 1] and FH_list[i] > FH_list[i + 1] and steps[i] > 5:
        max_x_FH.append(steps[i])
        max_y_FH.append(FH_list[i])

for i in range(1, len(FC_list) - 1):
    if FC_list[i] < FC_list[i - 1] and FC_list[i] < FC_list[i + 1]:
        min_x_FC.append(steps[i])
        min_y_FC.append(FC_list[i])
    if FH_list[i] < FH_list[i - 1] and FH_list[i] < FH_list[i + 1]:
        min_x_FH.append(steps[i])
        min_y_FH.append(FH_list[i])

print("Ширина головної пелюстки в площині H = " + str(round(SGP_FH, 2)) + '\u00b0')

len_value = [len(max_x_FC), len(max_x_FH), len(min_x_FC), len(min_x_FH)]
print("Табл. 2 - Аналіз ДС Рупора в площині H")
print("-----------------------------------------")
print("| № |θminC|θminH|θmaxC|FC(θ)|θmaxH|FH(θ)|")
for i in range(0, max(len_value)):
    v1 = f" {i + 1:.0f}"
    if len(v1) < len(str(max(len_value))) + 2: v1 += " "
    v2 = f"{min_x_FC[i]:.2f}" if i < len(min_x_FC) else '  -  '
    if len(v2) < 5: v2 += "0"
    v3 = f"{min_x_FH[i]:.2f}" if i < len(min_x_FH) else '  -  '
    if len(v3) < 5: v3 += "0"
    v4 = f"{max_x_FC[i]:.2f}" if i < len(max_x_FC) else '  -  '
    if len(v4) < 5: v4 += "0"
    h4 = f"{max_y_FC[i]:.3f}" if i < len(max_x_FC) else '  -  '
    if len(h4) < 5: h4 += "0"
    v5 = f"{max_x_FH[i]:.2f}" if i < len(max_x_FH) else '  -  '
    if len(v5) < 5: v5 += "0"
    h5 = f"{max_y_FH[i]:.3f}" if i < len(max_x_FH) else '  -  '
    if len(h5) < 5: h5 += "0"
    print(f"|{v1}|{v2}|{v3}|{v4}|{h4}|{v5}|{h5}|")
print("-----------------------------------------")

fig2, ax2 = plt.subplots(figsize=(20 / 2.54, 12 / 2.54))
ax2.plot(steps, F1_list, linewidth=0.7, label=r"$ F_{1h}(\theta) $ (Element)")
ax2.plot(steps, FC_list, linewidth=0.7, label=r"$ F_{CH}(\theta) $ (Aperture)")
ax2.plot(steps, FH_list, linewidth=0.7, label=r"$ F_{H}(\theta) $ (Total)")
ax2.plot(SGP_FH / 2, fS_FH, 'ro', markersize=4, label="ШГП в площині H")

ax2.annotate(f'({fS_FH:.3f}, {SGP_FH / 2:.2f}\u00b0)',
             xy=(SGP_FH / 2, fS_FH),
             xytext=((SGP_FH / 2) + 3, fS_FH + 0.05),
             arrowprops=dict(arrowstyle='->', color='black'), fontsize=8)

ax2.plot([0, SGP_FH / 2], [fS_FH, fS_FH], 'r--', linewidth=0.5)
ax2.plot([SGP_FH / 2, SGP_FH / 2], [0, fS_FH], 'r--', linewidth=0.5)

ax2.plot(max_x_FC, max_y_FC, "o", markersize=4, color="black", label=r"$ \theta_{max} $ FCH")
ax2.plot(max_x_FH, max_y_FH, "o", markersize=4, color="grey", label=r"$ \theta_{max} $ FH")
ax2.plot(min_x_FH, min_y_FH, "o", markersize=4, color="blue", label=r"$ \theta_{min} $ FH")

ax2.set_xlabel('θ' + '\u00b0', fontsize=10)
ax2.set_ylabel(r'|Fh($\theta$)|', fontsize=10)
plt.xticks(numpy.arange(0, 100, 2), fontsize=7)
plt.yticks(numpy.arange(0, 1.2, 0.1), fontsize=7)
plt.ylim(-0.01, 1.01)
plt.xlim(0, 90.5)
plt.title("Діаграма спрямованості в площині H")
plt.legend(loc="upper right", fontsize=7)
plt.grid(which='both', linestyle='--', linewidth=0.2, color='gray')
fig2.savefig("DS_Plane_H.jpg", dpi=600)

plt.show()