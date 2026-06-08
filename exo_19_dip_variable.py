import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import optuna
from scipy.optimize import differential_evolution
from scipy.spatial.distance import cdist
from scipy.signal import savgol_filter
import time

optuna.logging.set_verbosity(optuna.logging.WARNING)

# ==============================================================================
# --- 1. PARAMETROS ANTROPOMETRICOS (longitudes de falanges reales) ---
# ==============================================================================
FP_REAL = 0.049   # Falange proximal (m)
FM_REAL = 0.026   # Falange medial (m)
FD_REAL = 0.024   # Falange distal (m)

# ==============================================================================
# --- 2. CARGA Y CORRECCION DE DATOS MOCAP ---
# ==============================================================================
# Prioridad: datos de pinch (mayor rango DIP) > datos originales
# ==============================================================================
try:
    print(">> Intentando cargar 'mocap_pinch_pizzolato.csv' (datos de pinch)...")
    datos_mocap = pd.read_csv("mocap_pinch_pizzolato.csv")
    MOCAP_SOURCE = "mocap_pinch_pizzolato.csv"
    print(">> Exito: datos de pinch cargados.")
except FileNotFoundError:
    print(">> 'mocap_pinch_pizzolato.csv' no encontrado. Fallback a datos originales...")
    try:
        datos_mocap = pd.read_csv("mocap_indice_120pts.csv")
        MOCAP_SOURCE = "mocap_indice_120pts.csv"
        print(">> Fallback: 'mocap_indice_120pts.csv' cargado.")
    except FileNotFoundError:
        print("\n>> ERROR CRITICO: No se encontro ningun archivo MOCAP.")
        raise

# Procesamiento de los datos MOCAP
mcp_raw = datos_mocap['Theta_MCP'].values
pip_raw = datos_mocap['Theta_PIP'].values
dip_raw = datos_mocap['Theta_DIP'].values

if MOCAP_SOURCE == "mocap_indice_120pts.csv":
    # Los datos originales van de flexion maxima a extension; invertir
    mcp_raw = mcp_raw[::-1]
    pip_raw = pip_raw[::-1]
    dip_raw = dip_raw[::-1]
else:
    # Para datos de pinch: la secuencia es apertura -> cierre -> apertura (ciclo).
    # Extraer solo la fase de cierre: desde el inicio hasta el pico de flexion PIP.
    peak_idx = np.argmax(pip_raw)
    if peak_idx > 0 and peak_idx < len(pip_raw) - 1:
        # Tomar desde un punto cercano al inicio de movimiento hasta el pico
        # Buscar el inicio real del movimiento (primer cambio significativo)
        pip_diff = np.diff(pip_raw[:peak_idx+1])
        start_idx = 0
        for i in range(len(pip_diff)):
            if pip_diff[i] > 0.5:  # mas de 0.5 grados de cambio
                start_idx = i
                break
        mcp_raw = mcp_raw[start_idx:peak_idx+1]
        pip_raw = pip_raw[start_idx:peak_idx+1]
        dip_raw = dip_raw[start_idx:peak_idx+1]
        print(f"   Fase de cierre extraida: indices {start_idx} a {peak_idx} ({len(mcp_raw)} puntos)")

    # Resamplear a un numero fijo de puntos para consistencia con theta_input
    from scipy.interpolate import interp1d
    n_original = len(mcp_raw)
    N_TARGET = 120
    if n_original != N_TARGET and n_original > 3:
        t_orig = np.linspace(0, 1, n_original)
        t_new  = np.linspace(0, 1, N_TARGET)
        mcp_raw = interp1d(t_orig, mcp_raw, kind='linear')(t_new)
        pip_raw = interp1d(t_orig, pip_raw, kind='linear')(t_new)
        dip_raw = interp1d(t_orig, dip_raw, kind='linear')(t_new)
        print(f"   Resampleado de {n_original} a {N_TARGET} puntos.")

N_PUNTOS = len(mcp_raw)

# Recortar DIP a 0 deg minimo (no puede existir hiperextension en el exo)
dip_raw = np.clip(dip_raw, 0.0, None)

# Suavizado Savitzky-Golay para eliminar oscilaciones del sensor
WIN = 15  # ventana impar
if N_PUNTOS >= WIN:
    mcp_smooth = savgol_filter(np.deg2rad(mcp_raw), window_length=WIN, polyorder=2)
    pip_smooth = savgol_filter(np.deg2rad(pip_raw), window_length=WIN, polyorder=2)
    dip_smooth = savgol_filter(np.deg2rad(dip_raw), window_length=WIN, polyorder=2)
    dip_smooth = np.clip(dip_smooth, 0.0, None)
else:
    mcp_smooth = np.deg2rad(mcp_raw)
    pip_smooth = np.deg2rad(pip_raw)
    dip_smooth = np.deg2rad(dip_raw)

# Cinematica directa del dedo (cadena de cuerpos rigidos desde la MCF)
seg_prox = mcp_smooth
seg_med  = mcp_smooth + pip_smooth
seg_dist = mcp_smooth + pip_smooth + dip_smooth

pxIFP_mocap = FP_REAL * np.cos(seg_prox)
pyIFP_mocap = FP_REAL * np.sin(seg_prox)
pxIFD_mocap = pxIFP_mocap + FM_REAL * np.cos(seg_med)
pyIFD_mocap = pyIFP_mocap + FM_REAL * np.sin(seg_med)
pxPF_mocap  = pxIFD_mocap + FD_REAL * np.cos(seg_dist)
pyPF_mocap  = pyIFD_mocap + FD_REAL * np.sin(seg_dist)

mocap_pts = {
    'ifp': np.column_stack((pxIFP_mocap, pyIFP_mocap)),
    'ifd': np.column_stack((pxIFD_mocap, pyIFD_mocap)),
    'tip': np.column_stack((pxPF_mocap,  pyPF_mocap))
}

# Rango de entrada de la manivela principal (0 -> 85 deg)
theta_input = np.linspace(0, np.deg2rad(85), N_PUNTOS)

print(f">> MOCAP cargado ({MOCAP_SOURCE}): {N_PUNTOS} puntos.")
print(f"   MCP: {np.rad2deg(mcp_smooth[0]):.1f} deg -> {np.rad2deg(mcp_smooth[-1]):.1f} deg")
print(f"   PIP: {np.rad2deg(pip_smooth[0]):.1f} deg -> {np.rad2deg(pip_smooth[-1]):.1f} deg")
print(f"   DIP: {np.rad2deg(dip_smooth[0]):.1f} deg -> {np.rad2deg(dip_smooth[-1]):.1f} deg")
print(f"   Rango DIP: {np.rad2deg(dip_smooth.max() - dip_smooth.min()):.1f} deg")

# ==============================================================================
# --- 3. FUNCIONES DE EVALUACION ---
# ==============================================================================
def chamfer_distance(curve_target, curve_sim):
    """Distancia de Chamfer bidireccional (metrica de forma)."""
    dists = cdist(curve_target, curve_sim)
    return np.mean(np.min(dists, axis=1)) + np.mean(np.min(dists, axis=0))

def optimal_rigid_transform(target, sim):
    """Transformacion rigida optima (rotacion + traslacion) entre dos nubes de puntos."""
    c_t = np.mean(target, axis=0)
    c_s = np.mean(sim,    axis=0)
    H   = (sim - c_s).T @ (target - c_t)
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[1, :] *= -1
        R = Vt.T @ U.T
    t = c_t - R @ c_s
    return R, t

def apply_transform(points, R, t):
    return (R @ points.T).T + t

def monotonicity_penalty(curve):
    """
    Penaliza inversiones de direccion en la trayectoria simulada.
    Devuelve un valor >= 0 (0 = curva perfectamente monotona en avance).
    """
    diffs = np.diff(curve, axis=0)
    arc_increments = np.linalg.norm(diffs, axis=1)
    total_arc = np.sum(arc_increments)

    if total_arc < 1e-9:
        return 0.0
    mean_dir = np.sum(diffs, axis=0) / (total_arc + 1e-12)
    mean_dir /= (np.linalg.norm(mean_dir) + 1e-12)
    projections = diffs @ mean_dir
    reversals = np.sum(np.clip(-projections, 0, None))
    return reversals

# ==============================================================================
# --- 4. SOLUCIONES MATEMATICAS ---
# ==============================================================================
def sol_5_barras(r1, r2, r3, r4, r5, theta1, theta2):
    den = r4*np.cos(theta2) - r1*np.cos(theta1) + 2*r3
    if np.abs(den) < 1e-4:
        return None
    e = (r1*np.sin(theta1) - r4*np.sin(theta2)) / den
    f = (2*(r1*r3*np.cos(theta1) + r3*r4*np.cos(theta2))
         - r1**2 + r2**2 + r4**2 - r5**2) / (2*den)
    d_ = e**2 + 1
    g  = 2*(e*f - e*r1*np.cos(theta1) + e*r3 - r1*np.sin(theta1))
    h  = (f**2 - 2*f*(r1*np.cos(theta1) - r3)
          - 2*r1*r3*np.cos(theta1) + r1**2 + r3**2 - r2**2)
    disc = g**2 - 4*d_*h
    if disc < 0:
        return None
    py = (-g + np.sqrt(disc)) / (2*d_)
    px = e*py + f
    return px, py

def solve_four_bar(a, b, c, d, theta2, theta1):
    k1 = a*np.cos(theta2) + d*np.cos(theta1)
    k2 = a*np.sin(theta2) + d*np.sin(theta1)
    k3 = k1**2 + k2**2 + c**2 - b**2
    A1 = -2*k1*c - k3
    B1 =  4*k2*c
    C1 =  2*k1*c - k3
    disc = B1**2 - 4*A1*C1
    if disc < 0:
        return None
    return 2*np.arctan((-B1 - np.sqrt(disc)) / (2*A1))

# ==============================================================================
# --- 5. MODELO CINEMATICO CON DIP VARIABLE (POLINOMIO CUADRATICO) ---
# ==============================================================================
# Cambio clave respecto a exo_18.py:
#   - Se elimina theta_aux_fd (1 parametro constante, indice 14)
#   - Se agregan 3 coeficientes polinomiales: coeff_fd_0, coeff_fd_1, coeff_fd_2
#   - theta_aux_fd(step) = c0 + c1*delta_pip + c2*delta_pip^2
#     donde delta_pip = theta_fm - theta_fp (flexion relativa PIP)
#   - Total: 19 parametros (17 - 1 + 3 = 19)
# ==============================================================================
def run_kinematics(p, th_input):
    (Bancada1, Bancada2, Link1, Link2, Link3, Link4, Link5,
     Link6, Link7, Link8, Link10, hsp, dsp,
     theta_aux_fm, coeff_fd_0, coeff_fd_1, coeff_fd_2,
     gear_ratio, theta_offset) = p

    # Validaciones basicas
    if gear_ratio <= 0:
        return None
    if min(p[:11]) <= 0.005:
        return None
    if hsp <= 0 or dsp <= 0:
        return None
    if FP_REAL - 2.0 * dsp <= 0.001:
        return None

    # Pre-calculos fijos
    r3_val   = Bancada1 / 2.0
    theta14B = np.pi / 2.0
    c1       = np.sqrt(hsp**2 + dsp**2)

    rs2          = np.sqrt(hsp**2 + (FP_REAL - dsp)**2)
    theta_aux_s2 = np.arctan2(hsp, FP_REAL - dsp)

    PXifp, PYifp = [], []
    PXifd, PYifd = [], []
    PXtip, PYtip = [], []

    prev_theta_fm = None
    prev_theta_fd = None

    for th2 in th_input:
        th1 = (th2 / gear_ratio) + theta_offset

        # --- Mecanismo 1: 5 barras ---
        res5 = sol_5_barras(Link4, Link3, r3_val, Link1, Link2, th1, th2)
        if res5 is None:
            return None
        pxP, pyP = res5

        # --- Mecanismo 1: 4 barras ---
        theta4 = solve_four_bar(Link4, Link5, c1, Bancada2, th1, theta14B)
        if theta4 is None:
            return None

        # Angulo y posicion de la articulacion IFP
        theta_fp = theta4 + np.arctan2(hsp, dsp)
        px_ifp   = FP_REAL * np.cos(theta_fp) - Bancada2*np.cos(theta14B) - r3_val
        py_ifp   = FP_REAL * np.sin(theta_fp) - Bancada2*np.sin(theta14B)

        # Soportes de la falange proximal
        pxs1 = c1 * np.cos(theta4) - Bancada2*np.cos(theta14B) - r3_val
        pys1 = c1 * np.sin(theta4) - Bancada2*np.sin(theta14B)

        theta_ps2 = theta_fp - theta_aux_s2
        pxs2 = rs2 * np.cos(theta_ps2) - Bancada2*np.cos(theta14B) - r3_val
        pys2 = rs2 * np.sin(theta_ps2) - Bancada2*np.sin(theta14B)

        pxm4 = Link4*np.cos(th2) - r3_val
        pym4 = Link4*np.sin(th2)

        # --- Mecanismo 2: 5 barras (sistema de referencia secundario) ---
        theta_roll = np.arctan2(pym4 - pys1, pxm4 - pxs1)
        theta1m2   = np.arctan2(pys2 - pys1, pxs2 - pxs1) - theta_roll
        theta2m2   = np.arctan2(pyP  - pym4, pxP  - pxm4) - theta_roll

        res5_2 = sol_5_barras(FP_REAL - 2.0 * dsp, Link7, Link5/2.0, Link3, Link6,
                              theta1m2, theta2m2)
        if res5_2 is None:
            return None
        px_local, py_local = res5_2

        mag        = np.sqrt(px_local**2 + py_local**2)
        theta_loc  = np.arctan2(py_local, px_local)
        px_aux     = (pxs1 + pxm4) / 2.0
        py_aux     = (pys1 + pym4) / 2.0
        pxP2 = mag * np.cos(theta_loc + theta_roll) + px_aux
        pyP2 = mag * np.sin(theta_loc + theta_roll) + py_aux

        # --- Mecanismo 2: 4 barras ---
        theta1m42 = np.arctan2(pys2 - py_ifp, pxs2 - px_ifp)
        theta2m42 = np.arctan2(pyP2 - pys2,   pxP2 - pxs2)

        theta4m2 = solve_four_bar(Link7, Link8, Link10, c1, theta2m42, theta1m42)
        if theta4m2 is None:
            return None

        # Angulo de la falange medial
        theta_fm = theta4m2 + theta_aux_fm

        # ===============================================================
        # DIP VARIABLE: Acoplamiento polinomial
        # delta_pip = flexion relativa PIP (crece al cerrar el dedo)
        # theta_aux_fd = c0 + c1*delta_pip + c2*delta_pip^2
        # ===============================================================
        delta_pip = theta_fm - theta_fp
        theta_aux_fd = coeff_fd_0 + coeff_fd_1 * delta_pip + coeff_fd_2 * delta_pip**2
        theta_fd = theta_fm + theta_aux_fd

        # -----------------------------------------------------------------
        # RESTRICCION FISICA ANTI-GANCHO (falange medial):
        # theta_fm debe ser monotono durante el cierre.
        # -----------------------------------------------------------------
        if prev_theta_fm is not None:
            delta_fm = (theta_fm - prev_theta_fm + np.pi) % (2*np.pi) - np.pi
            if delta_fm < -np.deg2rad(0.5):
                return None
        prev_theta_fm = theta_fm

        # -----------------------------------------------------------------
        # RESTRICCION MONOTONIA DIP:
        # theta_fd debe ser monotono creciente (flexion progresiva).
        # Toleramos oscilaciones menores a 1 grado.
        # -----------------------------------------------------------------
        if prev_theta_fd is not None:
            delta_fd = (theta_fd - prev_theta_fd + np.pi) % (2*np.pi) - np.pi
            if delta_fd < -np.deg2rad(1.0):
                return None
        prev_theta_fd = theta_fd

        # Posiciones cartesianas
        px_ifd = FM_REAL * np.cos(theta_fm) + px_ifp
        py_ifd = FM_REAL * np.sin(theta_fm) + py_ifp
        px_tip = FD_REAL * np.cos(theta_fd) + px_ifd
        py_tip = FD_REAL * np.sin(theta_fd) + py_ifd

        # Rechazar valores fuera del espacio de trabajo
        if (np.abs(px_tip) > 0.3 or np.abs(py_tip) > 0.3
                or np.isnan(px_tip) or np.isnan(py_tip)):
            return None

        PXifp.append(px_ifp); PYifp.append(py_ifp)
        PXifd.append(px_ifd); PYifd.append(py_ifd)
        PXtip.append(px_tip); PYtip.append(py_tip)

    return {
        'ifp': np.column_stack((PXifp, PYifp)),
        'ifd': np.column_stack((PXifd, PYifd)),
        'tip': np.column_stack((PXtip, PYtip))
    }

# ==============================================================================
# --- 6. FUNCION OBJETIVO ---
# ==============================================================================
W_IFP = 0.25
W_IFD = 0.375
W_TIP = 0.375

W_MONO = 5.0

def fitness_function(p):
    sim_data = run_kinematics(p, theta_input)
    if sim_data is None:
        return 1000.0

    all_mocap = np.vstack([mocap_pts[k] for k in ('ifp', 'ifd', 'tip')])
    all_sim   = np.vstack([sim_data[k]  for k in ('ifp', 'ifd', 'tip')])

    R, t = optimal_rigid_transform(all_mocap, all_sim)

    aligned = {k: apply_transform(sim_data[k], R, t) for k in ('ifp', 'ifd', 'tip')}

    err_ifp = chamfer_distance(mocap_pts['ifp'], aligned['ifp'])
    err_ifd = chamfer_distance(mocap_pts['ifd'], aligned['ifd'])
    err_tip = chamfer_distance(mocap_pts['tip'], aligned['tip'])

    # Penalizacion de monotonicidad
    mono_ifd = monotonicity_penalty(aligned['ifd'])
    mono_tip = monotonicity_penalty(aligned['tip'])

    shape_error = W_IFP*err_ifp + W_IFD*err_ifd + W_TIP*err_tip
    mono_error  = W_MONO * (mono_ifd + mono_tip)

    return shape_error + mono_error

# ==============================================================================
# --- 7. BOUNDS (19 parametros) ---
# ==============================================================================
# Parametros 0-13: iguales a exo_18 (hasta theta_aux_fm)
# Parametros 14-16: coeff_fd_0, coeff_fd_1, coeff_fd_2 (reemplazan theta_aux_fd)
# Parametros 17-18: gear_ratio, theta_offset (iguales a exo_18 indices 15-16)
# ==============================================================================
bounds = [
    (0.015, 0.12), (0.015, 0.12),               # Bancada1, Bancada2
    (0.01,  0.12), (0.01,  0.12), (0.01, 0.12),  # Link1, Link2, Link3
    (0.01,  0.12), (0.01,  0.12),                 # Link4, Link5
    (0.01,  0.12), (0.01,  0.12), (0.01, 0.12),  # Link6, Link7, Link8
    (0.01,  0.12),                                 # Link10
    (-0.04, 0.08), (0.005, 0.023),                 # hsp, dsp
    (0.0,   np.pi),                                # theta_aux_fm
    (0.0,   np.pi),                                # coeff_fd_0 (offset base, similar a theta_aux_fd)
    (-2.0,  2.0),                                  # coeff_fd_1 (acoplamiento lineal)
    (-5.0,  5.0),                                  # coeff_fd_2 (acoplamiento cuadratico)
    (1.0,   8.0),                                  # gear_ratio
    (-np.pi, np.pi)                                # theta_offset
]

# ==============================================================================
# --- 8. OPTIMIZACION BAYESIANA (Optuna) ---
# ==============================================================================
def objective_optuna(trial):
    strategy      = trial.suggest_categorical('strategy',
                        ['best1exp', 'rand1bin', 'best1bin'])
    popsize       = trial.suggest_int('popsize', 15, 30)
    mut_min       = trial.suggest_float('mut_min', 0.5, 0.9)
    mut_max       = trial.suggest_float('mut_max', mut_min + 0.1, 1.5)
    recombination = trial.suggest_float('recombination', 0.5, 0.9)

    res = differential_evolution(
        fitness_function, bounds,
        strategy=strategy, popsize=popsize,
        mutation=(mut_min, mut_max), recombination=recombination,
        maxiter=20, tol=1e-3, polish=False,
        updating='immediate', workers=1
    )
    return res.fun

# ==============================================================================
# --- 9. EJECUCION PRINCIPAL ---
# ==============================================================================
if __name__ == '__main__':
    # =========================================================================
    # Validacion de retrocompatibilidad:
    # coeff_fd_0 = deg2rad(38.78), coeff_fd_1 = 0, coeff_fd_2 = 0
    # equivale al theta_aux_fd constante de exo_18.py
    # =========================================================================
    p_test = [0.018, 0.020, 0.035, 0.049, 0.025, 0.020, 0.025,
              0.055, 0.035, 0.052, 0.04601,
              0.017, 0.018,
              np.deg2rad(51.39),         # theta_aux_fm
              np.deg2rad(38.78),         # coeff_fd_0 (= antiguo theta_aux_fd)
              0.0,                        # coeff_fd_1 (sin acoplamiento lineal)
              0.0,                        # coeff_fd_2 (sin acoplamiento cuadratico)
              2.0,                        # gear_ratio
              np.deg2rad(109)]            # theta_offset

    print("\n>> VALIDACION DE RETROCOMPATIBILIDAD")
    print("   Parametros: coeff_fd = [38.78 deg, 0, 0] (equivalente a exo_18)")
    test_result = run_kinematics(p_test, theta_input)
    if test_result is not None:
        print("   RESULTADO: Cinematica OK - solucion valida")
        test_fitness = fitness_function(p_test)
        print(f"   Fitness: {test_fitness:.6f}")
        print(f"   Puntos generados: {len(test_result['tip'])}")

        # Verificar que la punta se mueve (no es un punto fijo)
        tip_range_x = np.ptp(test_result['tip'][:, 0]) * 1000
        tip_range_y = np.ptp(test_result['tip'][:, 1]) * 1000
        print(f"   Rango punta: X={tip_range_x:.2f} mm, Y={tip_range_y:.2f} mm")
    else:
        print("   RESULTADO: FALLO - la cinematica devolvio None")
        print("   NOTA: Esto puede ocurrir con datos de pinch que tienen mayor rango.")
        print("         Los parametros MATLAB fueron calibrados para los datos originales.")

    # =========================================================================
    # Prueba con coeficientes no triviales (acoplamiento activo)
    # =========================================================================
    print("\n>> PRUEBA CON ACOPLAMIENTO DIP-PIP ACTIVO")
    p_active = p_test.copy()
    p_active[14] = np.deg2rad(30.0)   # coeff_fd_0: offset base menor
    p_active[15] = 0.3                 # coeff_fd_1: acoplamiento lineal positivo
    p_active[16] = 0.1                 # coeff_fd_2: ligera curvatura

    test_active = run_kinematics(p_active, theta_input)
    if test_active is not None:
        print("   RESULTADO: Cinematica con acoplamiento activo - OK")
        fitness_active = fitness_function(p_active)
        print(f"   Fitness: {fitness_active:.6f}")
    else:
        print("   RESULTADO: Cinematica con acoplamiento activo - No converge")
        print("   (Normal: requiere optimizacion conjunta de todos los parametros)")

    # =========================================================================
    # Optimizacion completa (descomentariar para ejecutar)
    # =========================================================================
    print('\n====================================================')
    print('>> ETAPA 1: Busqueda de Estrategia con Optuna')
    print('====================================================')

    study = optuna.create_study(direction='minimize')
    t0 = time.time()
    study.optimize(objective_optuna, n_trials=8)
    best_hp = study.best_params
    print(f'\n>> Optuna finalizado en {time.time()-t0:.1f}s')
    print(f'   Mejor estrategia : {best_hp["strategy"]}')
    print(f'   Mejor popsize    : {best_hp["popsize"]}')

    print('\n====================================================')
    print('>> ETAPA 2: Optimizacion Cinematica Profunda (ED)')
    print('====================================================')

    resultado = differential_evolution(
        fitness_function, bounds,
        strategy=best_hp['strategy'],
        popsize=best_hp['popsize'],
        mutation=(best_hp['mut_min'], best_hp['mut_max']),
        recombination=best_hp['recombination'],
        maxiter=1000, tol=1e-5,
        polish=True, disp=True,
        updating='immediate', workers=1
    )

    p_opt    = resultado.x
    best_sim = run_kinematics(p_opt, theta_input)

    if best_sim is None:
        print("\n>> ADVERTENCIA: El resultado final no produce una cinematica valida.")
        print("   Considere ampliar los bounds o aumentar maxiter.")
        exit()

    # --- Transformacion de alineacion final ---
    all_mocap = np.vstack([mocap_pts[k] for k in ('ifp', 'ifd', 'tip')])
    all_sim   = np.vstack([best_sim[k]  for k in ('ifp', 'ifd', 'tip')])

    R_opt, t_opt = optimal_rigid_transform(all_mocap, all_sim)
    angulo_montaje = np.rad2deg(np.arctan2(R_opt[1, 0], R_opt[0, 0]))

    sim_aligned = {k: apply_transform(best_sim[k], R_opt, t_opt)
                   for k in ('ifp', 'ifd', 'tip')}

    err_ifp_mm = chamfer_distance(mocap_pts['ifp'], sim_aligned['ifp']) * 1000
    err_ifd_mm = chamfer_distance(mocap_pts['ifd'], sim_aligned['ifd']) * 1000
    err_tip_mm = chamfer_distance(mocap_pts['tip'], sim_aligned['tip']) * 1000
    error_global_mm = (err_ifp_mm + err_ifd_mm + err_tip_mm) / 3.0

    # --- Reporte de resultados ---
    print('\n====================================================')
    print('>> RESULTADOS CINEMATICOS (DIP VARIABLE)')
    print('====================================================')
    print(f"   Error Global : {error_global_mm:.3f} mm")
    print(f"   Error IFP    : {err_ifp_mm:.3f} mm")
    print(f"   Error IFD    : {err_ifd_mm:.3f} mm")
    print(f"   Error Punta  : {err_tip_mm:.3f} mm\n")

    nombres = [
        "Bancada1 (m)",            "Bancada2 (m)",
        "Link1 (m)",               "Link2 (m)",
        "Link3 (m)",               "Link4 (m)",
        "Link5 (m)",               "Link6 (m)",
        "Link7 (m)",               "Link8 (m)",
        "Link10 (m)",
        "hsp (m)",                 "dsp (m)",
        "Theta Aux FM (rad)",
        "Coeff FD 0 (rad)",        "Coeff FD 1",
        "Coeff FD 2",
        "Relacion de engranaje",   "Theta Offset (rad)"
    ]

    print('>> PARAMETROS DIMENSIONALES DEL MECANISMO (19 params)')
    for nombre, valor in zip(nombres, p_opt):
        print(f"   {nombre:30s}: {valor:.6f}")

    # Mostrar la funcion polinomial resultante
    c0, c1_coeff, c2_coeff = p_opt[14], p_opt[15], p_opt[16]
    print(f"\n>> FUNCION DE ACOPLAMIENTO DIP:")
    print(f"   theta_aux_fd = {np.rad2deg(c0):.2f} deg"
          f" + {c1_coeff:.4f}*delta_pip"
          f" + {c2_coeff:.4f}*delta_pip^2")
    print(f"   (delta_pip = theta_fm - theta_fp, en radianes)")

    print('\n>> PARAMETROS DE MONTAJE (transformacion rigida)')
    print(f'   Traslacion X : {t_opt[0]*1000:.2f} mm')
    print(f'   Traslacion Y : {t_opt[1]*1000:.2f} mm')
    print(f'   Rotacion Base: {angulo_montaje:.2f} deg')
    print('====================================================\n')

    # --- Grafica de resultados ---
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    ax = axes[0]
    lw = 2.5
    ax.plot(mocap_pts['ifp'][:, 0]*1000, mocap_pts['ifp'][:, 1]*1000,
            'r--', lw=lw, alpha=0.65, label='MOCAP IFP')
    ax.plot(mocap_pts['ifd'][:, 0]*1000, mocap_pts['ifd'][:, 1]*1000,
            'g--', lw=lw, alpha=0.65, label='MOCAP IFD')
    ax.plot(mocap_pts['tip'][:, 0]*1000, mocap_pts['tip'][:, 1]*1000,
            'b--', lw=lw, alpha=0.65, label='MOCAP Punta')

    ax.plot(sim_aligned['ifp'][:, 0]*1000, sim_aligned['ifp'][:, 1]*1000,
            'r-', lw=lw, label='EXO IFP')
    ax.plot(sim_aligned['ifd'][:, 0]*1000, sim_aligned['ifd'][:, 1]*1000,
            'g-', lw=lw, label='EXO IFD')
    ax.plot(sim_aligned['tip'][:, 0]*1000, sim_aligned['tip'][:, 1]*1000,
            'b-', lw=lw, label='EXO Punta')

    ax.set_aspect('equal')
    ax.grid(True, linestyle=':', alpha=0.7)
    ax.legend(loc='upper right', fontsize=9)
    ax.set_title(
        f'Biofidelidad con DIP Variable\nError Global: {error_global_mm:.3f} mm',
        fontsize=12, fontweight='bold')
    ax.set_xlabel('Eje X (mm)', fontsize=11)
    ax.set_ylabel('Eje Y (mm)', fontsize=11)

    # Grafica del polinomio de acoplamiento
    ax2 = axes[1]
    delta_range = np.linspace(-0.5, 2.0, 100)
    theta_aux_curve = c0 + c1_coeff * delta_range + c2_coeff * delta_range**2
    ax2.plot(np.rad2deg(delta_range), np.rad2deg(theta_aux_curve), 'k-', lw=2)
    ax2.axhline(y=38.78, color='r', linestyle='--', alpha=0.7,
                label='exo_18 constante (38.78 deg)')
    ax2.set_xlabel('delta_pip (deg)', fontsize=11)
    ax2.set_ylabel('theta_aux_fd (deg)', fontsize=11)
    ax2.set_title('Funcion de Acoplamiento DIP-PIP', fontsize=12, fontweight='bold')
    ax2.grid(True, linestyle=':', alpha=0.7)
    ax2.legend(fontsize=9)

    plt.tight_layout()
    plt.savefig('Resultados_DIP_Variable.png', dpi=150)
    plt.close()

    # --- Guardado ---
    np.savetxt("Parametros_DIP_Variable.txt", p_opt,
               header="19 parametros optimizados (DIP variable, polinomio cuadratico)")

    if len(mocap_pts['ifp']) == len(sim_aligned['ifp']):
        pd.DataFrame({
            'x_ifp_mocap': mocap_pts['ifp'][:, 0],
            'y_ifp_mocap': mocap_pts['ifp'][:, 1],
            'x_ifp_exo':   sim_aligned['ifp'][:, 0],
            'y_ifp_exo':   sim_aligned['ifp'][:, 1],
            'x_ifd_mocap': mocap_pts['ifd'][:, 0],
            'y_ifd_mocap': mocap_pts['ifd'][:, 1],
            'x_ifd_exo':   sim_aligned['ifd'][:, 0],
            'y_ifd_exo':   sim_aligned['ifd'][:, 1],
            'x_tip_mocap': mocap_pts['tip'][:, 0],
            'y_tip_mocap': mocap_pts['tip'][:, 1],
            'x_tip_exo':   sim_aligned['tip'][:, 0],
            'y_tip_exo':   sim_aligned['tip'][:, 1],
        }).to_csv('Resultados_DIP_Variable.csv', index=False)
        print(">> Archivos guardados: Parametros_DIP_Variable.txt, "
              "Resultados_DIP_Variable.csv, Resultados_DIP_Variable.png")
