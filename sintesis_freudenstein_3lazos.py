#!/usr/bin/env python3
"""
Sintesis de Freudenstein de 3 lazos para exoesqueleto de dedo indice.
Extrae datos de agarre fino de la base de datos Pizzolato y diseña 3 mecanismos
de 4 barras en cascada: Motor->MCP, MCP->PIP, PIP->DIP.

Autor: Optimizacion Exo2
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
import os

# =============================================================================
# SECCION 1: Extraccion de datos de agarre fino de la base de datos Pizzolato
# =============================================================================

print("=" * 70)
print("SINTESIS DE FREUDENSTEIN - 3 LAZOS PARA EXOESQUELETO DE DEDO")
print("=" * 70)
print("\n[1] Extrayendo datos de agarre fino de la base de datos Pizzolato...")

# Ruta al archivo de datos crudos
data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '..', 'Base-de-Datos-Pizzoleto', 'data',
                         'Raw data', 'DB1', 'S1_A1_E2_glove.txt')

if not os.path.exists(data_path):
    data_path = '/projects/sandbox/Base-de-Datos-Pizzoleto/data/Raw data/DB1/S1_A1_E2_glove.txt'

data = np.loadtxt(data_path)
print(f"    Archivo cargado: {data.shape[0]} muestras x {data.shape[1]} columnas")

# Columnas del dedo indice (0-indexed): 4=MCP, 5=PIP, 6=DIP
mcp_raw = data[:, 4]
pip_raw = data[:, 5]
dip_raw = data[:, 6]

# Rangos globales de los sensores
mcp_min, mcp_max = 55.0, 181.0
pip_min, pip_max = 14.0, 178.0
dip_min, dip_max = 60.1, 237.0

# Rangos anatomicos maximos (grados)
mcp_range_max = 90.0
pip_range_max = 100.0
dip_range_max = 70.0

# Calibracion lineal: sensor -> grados
mcp_deg = (mcp_raw - mcp_min) / (mcp_max - mcp_min) * mcp_range_max
pip_deg = (pip_raw - pip_min) / (pip_max - pip_min) * pip_range_max
dip_deg = (dip_raw - dip_min) / (dip_max - dip_min) * dip_range_max

# Detectar segmento de agarre fino (tip-pinch): region con alta flexion DIP
seg_start = 42300
seg_end = 42800
seg_mcp = mcp_deg[seg_start:seg_end]
seg_pip = pip_deg[seg_start:seg_end]
seg_dip = dip_deg[seg_start:seg_end]

# Encontrar la fase ascendente (reposo a pico de DIP)
peak_idx = np.argmax(seg_dip)

# Extraer fase ascendente
asc_mcp = seg_mcp[:peak_idx + 1]
asc_pip = seg_pip[:peak_idx + 1]
asc_dip = seg_dip[:peak_idx + 1]

print(f"    Segmento de agarre fino: muestras {seg_start}-{seg_end}")
print(f"    Fase ascendente: {len(asc_mcp)} muestras")

# Suavizar con filtro uniforme
smooth_size = 7
asc_mcp_s = uniform_filter1d(asc_mcp, size=smooth_size)
asc_pip_s = uniform_filter1d(asc_pip, size=smooth_size)
asc_dip_s = uniform_filter1d(asc_dip, size=smooth_size)

# Usar la porcion activa del movimiento (desde donde DIP empieza a moverse)
# Esto elimina la zona muerta inicial y da perfiles mas suaves para la sintesis
start_offset = 100  # Inicio de la porcion activa
asc_mcp_active = asc_mcp_s[start_offset:]
asc_pip_active = asc_pip_s[start_offset:]
asc_dip_active = asc_dip_s[start_offset:]

# Downsample a 120 puntos
n_target = 120
indices = np.linspace(0, len(asc_mcp_active) - 1, n_target, dtype=int)
mcp_120 = asc_mcp_active[indices] - asc_mcp_active[indices][0]
pip_120 = asc_pip_active[indices] - asc_pip_active[indices][0]
dip_120 = asc_dip_active[indices] - asc_dip_active[indices][0]

# Asegurar no-negativos y monotonia
mcp_120 = np.maximum(mcp_120, 0)
pip_120 = np.maximum(pip_120, 0)
dip_120 = np.maximum(dip_120, 0)
mcp_120 = np.maximum.accumulate(mcp_120)
pip_120 = np.maximum.accumulate(pip_120)
dip_120 = np.maximum.accumulate(dip_120)

# Rangos totales
range_mcp = mcp_120[-1]
range_pip = pip_120[-1]
range_dip = dip_120[-1]

print(f"    Rango MCP: 0 a {range_mcp:.1f} grados")
print(f"    Rango PIP: 0 a {range_pip:.1f} grados")
print(f"    Rango DIP: 0 a {range_dip:.1f} grados")

# Guardar referencia con los rangos completos
df_ref = pd.DataFrame({
    'Theta_MCP': mcp_120,
    'Theta_PIP': pip_120,
    'Theta_DIP': dip_120
})
df_ref.to_csv('agarre_fino_referencia.csv', index=False)
print("    Archivo guardado: agarre_fino_referencia.csv")


# =============================================================================
# SECCION 2: Implementacion de sintesis de Freudenstein de 3 puntos
# =============================================================================

print("\n[2] Implementando sintesis de Freudenstein de 3 puntos...")


def freudenstein_synthesis(phi_pts, psi_pts, d):
    """
    Sintesis de Freudenstein con 3 puntos de precision.

    Ecuacion de Freudenstein:
        K1*cos(psi) - K2*cos(phi) + K3 = cos(phi - psi)

    donde K1 = d/a, K2 = d/c, K3 = (a^2 - b^2 + c^2 + d^2)/(2*a*c)

    Parametros:
        phi_pts: array de 3 angulos de entrada (radianes)
        psi_pts: array de 3 angulos de salida (radianes)
        d: longitud del eslabon fijo (mm)

    Retorna:
        dict con a, b, c, d, K1, K2, K3, o None si no es valido
    """
    A = np.zeros((3, 3))
    B = np.zeros(3)

    for i in range(3):
        phi = phi_pts[i]
        psi = psi_pts[i]
        A[i, 0] = np.cos(psi)
        A[i, 1] = -np.cos(phi)
        A[i, 2] = 1.0
        B[i] = np.cos(phi - psi)

    try:
        det = np.linalg.det(A)
        if abs(det) < 1e-10:
            return None
        K = np.linalg.solve(A, B)
    except np.linalg.LinAlgError:
        return None

    K1, K2, K3 = K[0], K[1], K[2]

    if abs(K1) < 1e-10 or abs(K2) < 1e-10:
        return None
    a = d / K1
    c = d / K2

    if a <= 0 or c <= 0:
        return None

    b_sq = a**2 + c**2 + d**2 - 2.0 * a * c * K3
    if b_sq <= 0:
        return None

    b = np.sqrt(b_sq)

    return {
        'a': a, 'b': b, 'c': c, 'd': d,
        'K1': K1, 'K2': K2, 'K3': K3
    }


def simulate_four_bar(K1, K2, K3, phi_array, psi_start):
    """
    Simula el mecanismo de 4 barras usando sustitucion de medio-angulo.
    Sigue la rama correcta de forma continua desde psi_start.

    Ecuacion: K1*cos(psi) - K2*cos(phi) + K3 - cos(phi-psi) = 0
    Con t = tan(psi/2): cuadratica en t.
    """
    psi_out = np.full_like(phi_array, np.nan)
    psi_current = psi_start

    for i, phi in enumerate(phi_array):
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)

        P = K1 - cos_phi
        Q = -sin_phi
        R = K2 * cos_phi - K3

        # (P+R)*t^2 - 2*Q*t - (P-R) = 0
        a_coef = P + R
        b_coef = -2.0 * Q
        c_coef = -(P - R)

        if abs(a_coef) < 1e-15:
            if abs(b_coef) < 1e-15:
                continue
            t = -c_coef / b_coef
            psi_sol = 2.0 * np.arctan(t)
            psi_out[i] = psi_sol
            psi_current = psi_sol
            continue

        disc = b_coef**2 - 4.0 * a_coef * c_coef
        if disc < 0:
            continue

        sqrt_disc = np.sqrt(disc)
        t1 = (-b_coef + sqrt_disc) / (2.0 * a_coef)
        t2 = (-b_coef - sqrt_disc) / (2.0 * a_coef)

        psi1 = 2.0 * np.arctan(t1)
        psi2 = 2.0 * np.arctan(t2)

        # Elegir la solucion mas cercana (continuidad de rama)
        d1 = min(abs(psi1 - psi_current),
                 abs(psi1 + 2*np.pi - psi_current),
                 abs(psi1 - 2*np.pi - psi_current))
        d2 = min(abs(psi2 - psi_current),
                 abs(psi2 + 2*np.pi - psi_current),
                 abs(psi2 - 2*np.pi - psi_current))

        if d1 <= d2:
            psi_out[i] = psi1
            psi_current = psi1
        else:
            psi_out[i] = psi2
            psi_current = psi2

    return psi_out


def check_grashof(a, b, c, d):
    """Verifica condicion de Grashof: s + l <= p + q"""
    links = sorted([a, b, c, d])
    return (links[0] + links[3]) <= (links[1] + links[2])


def compute_rms_error(input_deg, output_deg, K1, K2, K3, phi0, psi0):
    """
    Calcula el error RMS para un mecanismo de 4 barras.
    """
    phi_full = np.radians(input_deg + phi0)
    psi_start = np.radians(psi0)

    psi_sim = simulate_four_bar(K1, K2, K3, phi_full, psi_start)
    psi_sim_deg = np.degrees(psi_sim) - psi0

    valid = ~np.isnan(psi_sim_deg)
    n_valid = np.sum(valid)

    if n_valid < len(output_deg) * 0.5:
        return np.inf, psi_sim_deg, n_valid / len(output_deg)

    error = psi_sim_deg[valid] - output_deg[valid]
    rms = np.sqrt(np.mean(error**2))
    coverage = n_valid / len(output_deg)

    return rms, psi_sim_deg, coverage


# =============================================================================
# SECCION 3: Definicion de lazos
# =============================================================================

print("\n[3] Definiendo lazos del mecanismo...")

# Rango del motor: 0 a 125 grados
motor_range = 125.0

# Para un mecanismo cascadeado:
# Loop 1: Motor -> MCP (reduccion de angulo)
# Loop 2: MCP -> PIP (amplificacion de angulo)
# Loop 3: PIP -> DIP (relacion cercana a 1:1 pero con reduccion)
#
# Cada lazo se diseña con relacion I/O lineal (distribucion uniforme).
# El rango de DIP se ajusta para garantizar solucion viable del mecanismo
# mientras se mantiene > 50 grados de flexion (requisito del agarre fino).

# Para Loop 3: el rango DIP debe ser menor que PIP para evitar
# la degeneracion de la ecuacion de Freudenstein (ratio ~1:1 no tiene solucion).
# Usamos 60 grados para DIP (> 50 requerido, < 64.6 de PIP).
range_dip_design = min(range_dip, range_pip - 5.0)
range_dip_design = max(range_dip_design, 55.0)  # Minimo 55 grados

print(f"    Lazo 1: Motor (0->{motor_range:.0f} deg) -> MCP (0->{range_mcp:.1f} deg)")
print(f"    Lazo 2: MCP (0->{range_mcp:.1f} deg) -> PIP (0->{range_pip:.1f} deg)")
print(f"    Lazo 3: PIP (0->{range_pip:.1f} deg) -> DIP (0->{range_dip_design:.1f} deg)")

# Trayectorias lineales para cada lazo
loop1_input = np.linspace(0, motor_range, n_target)
loop1_output = np.linspace(0, range_mcp, n_target)

loop2_input = np.linspace(0, range_mcp, n_target)
loop2_output = np.linspace(0, range_pip, n_target)

loop3_input = np.linspace(0, range_pip, n_target)
loop3_output = np.linspace(0, range_dip_design, n_target)

# Puntos de precision: inicio, medio, final
idx_prec = [0, n_target // 2, n_target - 1]

loops = [
    {'name': 'Lazo 1 (Motor->MCP)', 'input': loop1_input, 'output': loop1_output,
     'd_range': np.arange(15, 51, 5, dtype=float)},
    {'name': 'Lazo 2 (MCP->PIP)', 'input': loop2_input, 'output': loop2_output,
     'd_range': np.arange(15, 56, 5, dtype=float)},
    {'name': 'Lazo 3 (PIP->DIP)', 'input': loop3_input, 'output': loop3_output,
     'd_range': np.arange(10, 36, 5, dtype=float)},
]


# =============================================================================
# SECCION 4: Optimizacion por busqueda exhaustiva
# =============================================================================

print("\n[4] Buscando dimensiones optimas para cada lazo...")

LINK_MIN = 3.0
LINK_MAX = 70.0


def search_loop(loop_def, link_min=LINK_MIN, link_max=LINK_MAX,
                phi0_step=5, psi0_step=5, max_ratio=6.0):
    """
    Busqueda exhaustiva de phi0, psi0 y d para un lazo.
    Evalua candidatos por error RMS de simulacion completa.
    """
    input_deg = loop_def['input']
    output_deg = loop_def['output']
    d_range = loop_def['d_range']

    delta_phi = np.array([input_deg[idx_prec[0]], input_deg[idx_prec[1]],
                          input_deg[idx_prec[2]]])
    delta_psi = np.array([output_deg[idx_prec[0]], output_deg[idx_prec[1]],
                          output_deg[idx_prec[2]]])

    best_solution = None
    best_rms = np.inf

    phi0_values = np.arange(5, 176, phi0_step)
    psi0_values = np.arange(5, 176, psi0_step)

    for phi0 in phi0_values:
        for psi0 in psi0_values:
            phi_pts = np.radians(delta_phi + phi0)
            psi_pts = np.radians(delta_psi + psi0)

            for d in d_range:
                result = freudenstein_synthesis(phi_pts, psi_pts, d)
                if result is None:
                    continue

                a, b, c = result['a'], result['b'], result['c']
                links = [a, b, c, d]

                if any(l < link_min for l in links):
                    continue
                if any(l > link_max for l in links):
                    continue

                ratio = max(links) / min(links)
                if ratio > max_ratio:
                    continue

                # Evaluar error RMS
                rms, _, coverage = compute_rms_error(
                    input_deg, output_deg,
                    result['K1'], result['K2'], result['K3'],
                    phi0, psi0
                )

                if rms < best_rms and coverage > 0.8:
                    best_rms = rms
                    best_solution = {
                        'a': a, 'b': b, 'c': c, 'd': d,
                        'K1': result['K1'], 'K2': result['K2'],
                        'K3': result['K3'],
                        'phi0': phi0, 'psi0': psi0,
                        'ratio': ratio,
                        'grashof': check_grashof(a, b, c, d),
                        'rms': rms,
                        'coverage': coverage
                    }

    return best_solution


results = []

for i, loop in enumerate(loops):
    print(f"\n    Optimizando {loop['name']}...")

    sol = search_loop(loop)

    if sol is None:
        print(f"      No se encontro con restricciones normales. Relajando...")
        sol = search_loop(loop, link_min=2.0, link_max=80.0, max_ratio=8.0)

    if sol is None:
        print(f"      Segundo intento con grid fino...")
        sol = search_loop(loop, link_min=2.0, link_max=80.0,
                          phi0_step=3, psi0_step=3, max_ratio=10.0)

    if sol is not None:
        print(f"      SOLUCION ENCONTRADA:")
        print(f"        a={sol['a']:.2f} mm, b={sol['b']:.2f} mm, "
              f"c={sol['c']:.2f} mm, d={sol['d']:.2f} mm")
        print(f"        phi0={sol['phi0']:.0f} deg, psi0={sol['psi0']:.0f} deg")
        print(f"        Ratio max/min = {sol['ratio']:.2f}")
        print(f"        Grashof: {'Si' if sol['grashof'] else 'No'}")
        print(f"        Error RMS = {sol['rms']:.3f} deg")
    else:
        print(f"      ERROR: No se encontro solucion para {loop['name']}")

    results.append(sol)


# =============================================================================
# SECCION 5: Simulacion directa y verificacion de error
# =============================================================================

print("\n[5] Verificando soluciones mediante simulacion directa...")

errors_rms = []
sim_outputs = []

for i, (loop, sol) in enumerate(zip(loops, results)):
    if sol is None:
        errors_rms.append(np.nan)
        sim_outputs.append(None)
        print(f"    {loop['name']}: SIN SOLUCION - no se puede simular")
        continue

    K1, K2, K3 = sol['K1'], sol['K2'], sol['K3']
    phi0, psi0 = sol['phi0'], sol['psi0']
    input_deg = loop['input']
    output_deg = loop['output']

    rms, psi_sim_deg, coverage = compute_rms_error(
        input_deg, output_deg, K1, K2, K3, phi0, psi0
    )

    errors_rms.append(rms)
    sim_outputs.append(psi_sim_deg)

    print(f"    {loop['name']}: Error RMS = {rms:.2f} deg, "
          f"Cobertura = {coverage*100:.1f}%")


# =============================================================================
# SECCION 6: Salida de resultados
# =============================================================================

print("\n[6] Generando archivos de salida...")

with open('resultados_freudenstein.txt', 'w') as f:
    f.write("=" * 70 + "\n")
    f.write("RESULTADOS DE SINTESIS DE FREUDENSTEIN - 3 LAZOS\n")
    f.write("Exoesqueleto de dedo indice - Agarre fino (tip pinch)\n")
    f.write("=" * 70 + "\n")
    f.write("\n")
    f.write("Datos de referencia extraidos de: Base de datos Pizzolato\n")
    f.write(f"  Rango MCP: 0 a {range_mcp:.1f} grados\n")
    f.write(f"  Rango PIP: 0 a {range_pip:.1f} grados\n")
    f.write(f"  Rango DIP (diseno): 0 a {range_dip_design:.1f} grados\n")
    f.write(f"  Rango DIP (dato): 0 a {range_dip:.1f} grados\n")
    f.write(f"  Rango Motor: 0 a {motor_range:.0f} grados\n")
    f.write("\n")
    f.write("Nota: Cada lazo mapea linealmente su entrada a la salida.\n")
    f.write("El mecanismo cascadeado reproduce las amplitudes del agarre fino.\n")
    f.write("\n")

    for i, (loop, sol) in enumerate(zip(loops, results)):
        f.write("-" * 50 + "\n")
        f.write(f"{loop['name']}\n")
        f.write("-" * 50 + "\n")
        if sol is not None:
            f.write(f"  Eslabon fijo (d):   {sol['d']:.2f} mm\n")
            f.write(f"  Manivela (a):       {sol['a']:.2f} mm\n")
            f.write(f"  Acoplador (b):      {sol['b']:.2f} mm\n")
            f.write(f"  Balancin (c):       {sol['c']:.2f} mm\n")
            f.write(f"  Angulo inicial entrada (phi0): {sol['phi0']:.0f} deg\n")
            f.write(f"  Angulo inicial salida  (psi0): {sol['psi0']:.0f} deg\n")
            f.write(f"  Ratio max/min:      {sol['ratio']:.2f}\n")
            f.write(f"  Condicion Grashof:  {'Si' if sol['grashof'] else 'No'}\n")
            if not np.isnan(errors_rms[i]):
                f.write(f"  Error RMS:          {errors_rms[i]:.2f} deg\n")
            f.write("\n")
        else:
            f.write("  ERROR: No se encontro solucion viable\n")
            f.write("\n")

    f.write("=" * 70 + "\n")
    f.write("Dimensiones del dedo:\n")
    f.write("  Falange proximal (fp): 49 mm\n")
    f.write("  Falange medial (fm):   26 mm\n")
    f.write("  Falange distal (fd):   24 mm\n")
    f.write("=" * 70 + "\n")

print("    Archivo guardado: resultados_freudenstein.txt")

# Generar grafica
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
loop_labels = ['MCP (Lazo 1)', 'PIP (Lazo 2)', 'DIP (Lazo 3)']
input_labels = ['Motor (deg)', 'MCP (deg)', 'PIP (deg)']

for i, (loop, sol, ax) in enumerate(zip(loops, results, axes)):
    input_deg = loop['input']
    output_deg = loop['output']

    ax.plot(input_deg, output_deg, 'b-', linewidth=2, label='Deseado')

    if sim_outputs[i] is not None:
        valid = ~np.isnan(sim_outputs[i])
        if np.any(valid):
            rms_label = f'Freudenstein (RMS={errors_rms[i]:.2f} deg)'
            ax.plot(input_deg[valid], sim_outputs[i][valid], 'r--',
                    linewidth=1.5, label=rms_label)

    ax.set_xlabel(input_labels[i])
    ax.set_ylabel(f'{loop_labels[i]}')
    ax.set_title(f'{loop["name"]}')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    if sol is not None:
        for idx in idx_prec:
            ax.plot(input_deg[idx], output_deg[idx], 'go', markersize=10, zorder=5)

plt.tight_layout()
plt.savefig('sintesis_freudenstein_resultado.png', dpi=150, bbox_inches='tight')
plt.close()
print("    Archivo guardado: sintesis_freudenstein_resultado.png")


# =============================================================================
# SECCION 7: Resumen final
# =============================================================================

print("\n" + "=" * 70)
print("RESUMEN FINAL")
print("=" * 70)
all_ok = True
for i, (loop, sol) in enumerate(zip(loops, results)):
    if sol is not None and not np.isnan(errors_rms[i]):
        status = "OK" if errors_rms[i] < 5.0 else "ADVERTENCIA: Error alto"
        print(f"  {loop['name']}: {status} (RMS={errors_rms[i]:.2f} deg)")
        if errors_rms[i] >= 5.0:
            all_ok = False
    else:
        print(f"  {loop['name']}: FALLO")
        all_ok = False

if all_ok:
    print("\n  SINTESIS COMPLETADA EXITOSAMENTE")
else:
    print("\n  SINTESIS COMPLETADA CON ADVERTENCIAS")

print("=" * 70)
print("\nArchivos generados:")
print("  - agarre_fino_referencia.csv")
print("  - resultados_freudenstein.txt")
print("  - sintesis_freudenstein_resultado.png")
print()
