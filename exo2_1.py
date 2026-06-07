"""
exo2_1.py
=========

Optimizacion cinematica del exoesqueleto de dedo indice CON el tercer
mecanismo de 4 barras (4B#3) que da movimiento INDEPENDIENTE a la
falange distal (DIP). Sustituye al exo_18.py (que asumia una relacion
rigida theta_fd = theta_fm + theta_aux_fd).

CAMBIOS RESPECTO A exo_18.py
----------------------------
1. **La cinematica usada AHORA es la validada contra
   `CinematicaExoModificada.m`** (ver `diagrama_mecanismo.compute_geometry`).
   exo_18.py tenia una formulacion paralela que no era exactamente
   equivalente al MATLAB en posiciones absolutas; ese script se sustituye.
2. El bloque del DIP usa el 4B#3:
   bancada = falange medial Fm,
   manivela Lpc rigida a Fp con offset BETA1,
   balancin Lpd rigida a Fd con offset BETA2,
   acoplador Lac flotante.
3. Vector de parametros: pasa de 17 a 22 entradas. Se eliminan los que el
   modelo viejo usaba para el offset rigido y se agregan Lpc, Lpd, Lac,
   BETA1, BETA2 (mas un fallback rigido por si el 4B#3 no ensambla).
4. Bounds de los nuevos parametros tomados de Analisis_Modificacion_DIP.md.
5. La continuidad de alpha2_3 se garantiza con desenrollado entre pasos
   consecutivos del bucle (igual que en CinematicaExoModificada.m).

CONVENCION DE ANGULOS BETA1, BETA2 (importante)
-----------------------------------------------
Hay DOS convenciones para los angulos de montaje de los postes Lpc y Lpd:

  (A) PARAMETROS ANALITICOS - los que optimiza este script:
      BETA1, BETA2 son los offsets de la formulacion del lado palmar
      (rama "-" del arccos). Para los valores tipicos (BETA1 ~40 deg,
      BETA2 ~110 deg) la trayectoria DIP sale 32 -> 62 deg monotona. NO
      son los angulos que se miden con transportador en el dibujo dorsal.

  (B) ANGULOS GEOMETRICOS - los que se usan en CAD:
      BETA1_geom, BETA2_geom son los angulos fisicos medidos entre el eje
      de la falange y el poste, en la representacion DORSAL real. Para
      los mismos valores tipicos resultan BETA1_geom ~ -11.4 deg y
      BETA2_geom ~ 185.8 deg. Estos son los angulos para construir el CAD.

Ambas convenciones describen el MISMO mecanismo fisico (ver
Analisis_Modificacion_DIP.md, seccion 5.1). Las trayectorias son
identicas. Al final del script se reportan los valores _geom calculados
a partir de los analiticos optimizados, listos para el plano del CAD.
"""

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
# --- 1. PARAMETROS ANTROPOMETRICOS ---
# ==============================================================================
FP_REAL = 0.049   # Falange proximal (m)
FM_REAL = 0.026   # Falange medial (m)
FD_REAL = 0.024   # Falange distal (m)


# ==============================================================================
# --- 2. CARGA Y CORRECCION DE DATOS MOCAP ---
# ==============================================================================
try:
    print(">> Cargando base de datos MOCAP 'mocap_indice_120pts.csv'...")
    datos_mocap = pd.read_csv("mocap_indice_120pts.csv")

    mcp_raw = datos_mocap['Theta_MCP'].values[::-1]
    pip_raw = datos_mocap['Theta_PIP'].values[::-1]
    dip_raw = datos_mocap['Theta_DIP'].values[::-1]
    N_PUNTOS = len(mcp_raw)

    dip_raw = np.clip(dip_raw, 0.0, None)

    WIN = 15
    if N_PUNTOS >= WIN:
        mcp_smooth = savgol_filter(np.deg2rad(mcp_raw), window_length=WIN, polyorder=2)
        pip_smooth = savgol_filter(np.deg2rad(pip_raw), window_length=WIN, polyorder=2)
        dip_smooth = savgol_filter(np.deg2rad(dip_raw), window_length=WIN, polyorder=2)
        dip_smooth = np.clip(dip_smooth, 0.0, None)
    else:
        mcp_smooth = np.deg2rad(mcp_raw)
        pip_smooth = np.deg2rad(pip_raw)
        dip_smooth = np.deg2rad(dip_raw)
except FileNotFoundError:
    print("\n>> ERROR CRITICO: No se encontro 'mocap_indice_120pts.csv'.")
    raise

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

# THETA2 va de 0 a 132 deg, igual que CinematicaExoModificada.m
theta_input = np.linspace(0, np.deg2rad(132), N_PUNTOS)

print(f">> MOCAP cargado: {N_PUNTOS} puntos.")
print(f"   MCP: {np.rad2deg(mcp_smooth[0]):.1f} -> {np.rad2deg(mcp_smooth[-1]):.1f} deg")
print(f"   PIP: {np.rad2deg(pip_smooth[0]):.1f} -> {np.rad2deg(pip_smooth[-1]):.1f} deg")
print(f"   DIP: {np.rad2deg(dip_smooth[0]):.1f} -> {np.rad2deg(dip_smooth[-1]):.1f} deg")


# ==============================================================================
# --- 3. METRICAS DE EVALUACION ---
# ==============================================================================
def chamfer_distance(curve_target, curve_sim):
    dists = cdist(curve_target, curve_sim)
    return np.mean(np.min(dists, axis=1)) + np.mean(np.min(dists, axis=0))


def optimal_rigid_transform(target, sim):
    c_t = np.mean(target, axis=0)
    c_s = np.mean(sim, axis=0)
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
# --- 4. CINEMATICA VALIDADA (mismas ecuaciones que CinematicaExoModificada.m) ---
# ==============================================================================
def four_bar_open(a_, b_, c_, d_, th2_, th1_):
    """4 barras configuracion abierta (mismas ecuaciones que el MATLAB)."""
    k1 = d_*np.cos(th1_) + a_*np.cos(th2_)
    k2 = d_*np.sin(th1_) + a_*np.sin(th2_)
    k3 = k1**2 + k2**2 + c_**2 - b_**2
    A1 = -k3 - 2*k1*c_
    B1 = 4*k2*c_
    C1 = 2*k1*c_ - k3
    disc = B1**2 - 4*A1*C1
    if disc < 0:
        return None
    return 2*np.arctan((-B1 - np.sqrt(disc)) / (2*A1))


def run_kinematics(p, th_input):
    """Cinematica directa del exoesqueleto, IDENTICA a CinematicaExoModificada.m
    (verificada contra `diagrama_mecanismo.compute_geometry`).

    Vector p (22 entradas, todas en SI: metros / radianes / adimensional):
        [ 0]  Bancada1   (m)   distancia entre pivotes A y B
        [ 1]  Bancada2   (m)   distancia A - MCF (vertical)
        [ 2]  L1         (m)   manivela 5B#1 desde B
        [ 3]  L2         (m)   eslabon J2-P del 5B#1
        [ 4]  L3         (m)   eslabon M4-P del 5B#1
        [ 5]  L4         (m)   manivela 5B#1 desde A / manivela 4B#1
        [ 6]  L5         (m)   eslabon M4-S1 del 4B#1
        [ 7]  L6         (m)   eslabon P-P2 del 5B#2
        [ 8]  L7         (m)   manivela 4B#2 (S2-P2)
        [ 9]  L8         (m)   acoplador 4B#2 (P2-P3)
        [10]  c2         (m)   balancin 4B#2 (IFP-P3) -- en MATLAB se llama c2 = 46.01 mm
        [11]  hsp        (m)   altura del soporte sobre la falange
        [12]  dsp        (m)   distancia del soporte a la articulacion
        [13]  theta_aux_fm (rad)  offset entre c2 y Fm
        [14]  gear_ratio (-)    relacion de engranaje
        [15]  theta_offset (rad) offset angular del engranaje
        [16]  Lpc        (m)   manivela 4B#3 (poste sobre Fp)
        [17]  Lpd        (m)   balancin 4B#3 (poste sobre Fd)
        [18]  Lac        (m)   acoplador 4B#3
        [19]  BETA1      (rad) **angulo analitico** de montaje Lpc
        [20]  BETA2      (rad) **angulo analitico** de montaje Lpd
        [21]  theta_aux_fd_fb (rad) fallback rigido si 4B#3 no ensambla

    Devuelve dict con 'ifp', 'ifd', 'tip' (cada uno Nx2 en metros) o None
    si el set de parametros no es valido.
    """
    (Bancada1, Bancada2, L1, L2, L3, L4, L5,
     L6, L7, L8, c2_par, hsp, dsp,
     theta_aux_fm, gear_ratio, theta_offset,
     Lpc, Lpd, Lac, BETA1, BETA2,
     theta_aux_fd_fb) = p

    # Validaciones tempranas
    if gear_ratio <= 0:
        return None
    if min(Bancada1, Bancada2, L1, L2, L3, L4, L5, L6, L7, L8, c2_par,
           Lpc, Lpd, Lac) <= 1e-3:
        return None
    if hsp <= 0 or dsp <= 0:
        return None
    if FP_REAL - 2.0 * dsp <= 1e-3:
        return None

    # --- Constantes derivadas (igual que en MATLAB) ---
    r1 = L4; r2 = L3; r3 = Bancada1 / 2.0; r4 = L1; r5 = L2
    a  = L4; b  = L5; c_arm = np.sqrt(hsp**2 + dsp**2); d_arm = Bancada2
    theta14B = np.pi / 2.0
    r1m2 = FP_REAL - 2.0 * dsp
    r2m2 = L7; r3m2 = L5 / 2.0; r4m2 = L3; r5m2 = L6
    a2 = L7; b2 = L8; d2 = np.sqrt(hsp**2 + dsp**2)
    rs2          = np.sqrt(hsp**2 + (FP_REAL - dsp)**2)
    theta_aux_s2 = np.arctan2(hsp, FP_REAL - dsp)

    # Marco fijo (mismas coords absolutas que diagrama_mecanismo.py)
    A   = np.array([-r3, 0.0])
    Bp  = np.array([+r3, 0.0])
    MCF = np.array([-r3, -d_arm])

    PXifp, PYifp = [], []
    PXifd, PYifd = [], []
    PXtip, PYtip = [], []

    prev_theta_fm = None
    alpha2_3_prev = None

    for th2 in th_input:
        th1 = (th2 / gear_ratio) + theta_offset

        # --- Mecanismo 1: 5 barras ---
        M4 = A + r1 * np.array([np.cos(th1), np.sin(th1)])
        J2 = Bp + r4 * np.array([np.cos(th2), np.sin(th2)])
        den = r4*np.cos(th2) - r1*np.cos(th1) + 2*r3
        if abs(den) < 1e-6:
            return None
        e_ = (r1*np.sin(th1) - r4*np.sin(th2)) / den
        f_ = (2*(r1*r3*np.cos(th1) + r3*r4*np.cos(th2))
              - r1**2 + r2**2 + r4**2 - r5**2) / (2*den)
        daux = e_**2 + 1
        g_ = 2*(e_*f_ - e_*r1*np.cos(th1) + e_*r3 - r1*np.sin(th1))
        h_ = (f_**2 - 2*f_*(r1*np.cos(th1) - r3)
              - 2*r1*r3*np.cos(th1) + r1**2 + r3**2 - r2**2)
        disc5 = g_**2 - 4*daux*h_
        if disc5 < 0:
            return None
        pyP = (-g_ + np.sqrt(disc5)) / (2*daux)
        pxP = e_*pyP + f_
        Pp  = np.array([pxP, pyP])

        # --- Mecanismo 1: 4 barras ---
        th4a = four_bar_open(a, b, c_arm, d_arm, th1, theta14B)
        if th4a is None:
            return None
        TH4a_deg = np.rad2deg(th4a)
        if TH4a_deg < 0:
            TH4a_deg += 360
        th4a = np.deg2rad(TH4a_deg)
        theta_fp = np.deg2rad(TH4a_deg + np.rad2deg(np.arctan2(hsp, dsp)))
        IFP = MCF + FP_REAL * np.array([np.cos(theta_fp), np.sin(theta_fp)])
        S1  = MCF + c_arm * np.array([np.cos(th4a), np.sin(th4a)])
        thps2 = theta_fp - theta_aux_s2
        S2 = MCF + rs2 * np.array([np.cos(thps2), np.sin(thps2)])

        # --- Mecanismo 2: 5 barras ---
        throll = np.arctan2(M4[1] - S1[1], M4[0] - S1[0])

        def ang_local(pt_to, pt_from):
            ang = np.rad2deg(np.arctan2(pt_to[1] - pt_from[1], pt_to[0] - pt_from[0]))
            if ang < 0:
                ang += 360
            return np.deg2rad(ang) - throll

        th1m2 = ang_local(S2, S1)
        th2m2 = ang_local(Pp, M4)
        denm2 = r4m2*np.cos(th2m2) - r1m2*np.cos(th1m2) + 2*r3m2
        if abs(denm2) < 1e-6:
            return None
        em2 = (r1m2*np.sin(th1m2) - r4m2*np.sin(th2m2)) / denm2
        fm2 = (2*(r1m2*r3m2*np.cos(th1m2) + r3m2*r4m2*np.cos(th2m2))
               - r1m2**2 + r2m2**2 + r4m2**2 - r5m2**2) / (2*denm2)
        dauxm2 = em2**2 + 1
        gm2 = 2*(em2*fm2 - em2*r1m2*np.cos(th1m2) + em2*r3m2 - r1m2*np.sin(th1m2))
        hm2 = (fm2**2 - 2*fm2*(r1m2*np.cos(th1m2) - r3m2)
               - 2*r1m2*r3m2*np.cos(th1m2) + r1m2**2 + r3m2**2 - r2m2**2)
        discm2 = gm2**2 - 4*dauxm2*hm2
        if discm2 < 0:
            return None
        pyP2 = (-gm2 + np.sqrt(discm2)) / (2*dauxm2)
        pxP2 = em2*pyP2 + fm2
        p2_  = np.hypot(pxP2, pyP2)
        th2p2 = np.arctan2(pyP2, pxP2)
        AUX = (S1 + M4) / 2.0
        P2 = p2_ * np.array([np.cos(th2p2 + throll),
                             np.sin(th2p2 + throll)]) + AUX

        # --- Mecanismo 2: 4 barras (driver de Fm) ---
        th14B2 = np.arctan2(S2[1] - IFP[1], S2[0] - IFP[0])
        a24_deg = np.rad2deg(np.arctan2(P2[1] - S2[1], P2[0] - S2[0]))
        if a24_deg < 0:
            a24_deg += 360
        th24B2 = np.deg2rad(a24_deg)
        th4am2 = four_bar_open(a2, b2, c2_par, d2, th24B2, th14B2)
        if th4am2 is None:
            return None
        TH4am2_deg = np.rad2deg(th4am2)
        if TH4am2_deg < 0:
            TH4am2_deg += 360
        theta_fm = np.deg2rad(TH4am2_deg + np.rad2deg(theta_aux_fm))
        IFD = IFP + FM_REAL * np.array([np.cos(theta_fm), np.sin(theta_fm)])

        # ----- restriccion anti-gancho en la medial -----
        if prev_theta_fm is not None:
            delta_fm = (theta_fm - prev_theta_fm + np.pi) % (2*np.pi) - np.pi
            if delta_fm < -np.deg2rad(0.5):
                return None
        prev_theta_fm = theta_fm

        # --- Mecanismo 3: 4 barras (DIP) ---
        alpha1_3 = theta_fp + BETA1 - theta_fm
        Ax3 = Lpc * np.cos(alpha1_3)
        Ay3 = Lpc * np.sin(alpha1_3)
        Px3 = Ax3 - FM_REAL
        Py3 = Ay3
        R3  = np.hypot(Px3, Py3)
        K3  = (Px3**2 + Py3**2 + Lpd**2 - Lac**2) / (2.0 * Lpd)

        if abs(K3) > R3:
            theta_fd = theta_fm + theta_aux_fd_fb
            alpha2_3_prev = None
        else:
            phi3 = np.arctan2(Py3, Px3)
            alpha2_3 = phi3 - np.arccos(np.clip(K3 / R3, -1.0, 1.0))
            if alpha2_3_prev is not None:
                while alpha2_3 - alpha2_3_prev > np.pi:
                    alpha2_3 -= 2*np.pi
                while alpha2_3 - alpha2_3_prev < -np.pi:
                    alpha2_3 += 2*np.pi
            alpha2_3_prev = alpha2_3
            theta_fd = theta_fm + (alpha2_3 - BETA2)

        TIP = IFD + FD_REAL * np.array([np.cos(theta_fd), np.sin(theta_fd)])

        if (np.abs(TIP[0]) > 0.3 or np.abs(TIP[1]) > 0.3
                or np.isnan(TIP[0]) or np.isnan(TIP[1])):
            return None

        PXifp.append(IFP[0]); PYifp.append(IFP[1])
        PXifd.append(IFD[0]); PYifd.append(IFD[1])
        PXtip.append(TIP[0]); PYtip.append(TIP[1])

    return {
        'ifp': np.column_stack((PXifp, PYifp)),
        'ifd': np.column_stack((PXifd, PYifd)),
        'tip': np.column_stack((PXtip, PYtip))
    }


# ==============================================================================
# --- 5. FUNCION OBJETIVO ---
# ==============================================================================
W_IFP  = 0.25
W_IFD  = 0.375
W_TIP  = 0.375
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
    mono_ifd = monotonicity_penalty(aligned['ifd'])
    mono_tip = monotonicity_penalty(aligned['tip'])
    return (W_IFP*err_ifp + W_IFD*err_ifd + W_TIP*err_tip
            + W_MONO*(mono_ifd + mono_tip))


# ==============================================================================
# --- 6. BOUNDS (22 parametros) ---
# ==============================================================================
bounds = [
    (0.015, 0.12), (0.015, 0.12),                  # Bancada1, Bancada2
    (0.01, 0.12), (0.01, 0.12), (0.01, 0.12),      # L1, L2, L3
    (0.01, 0.12), (0.01, 0.12),                     # L4, L5
    (0.01, 0.12), (0.01, 0.12), (0.01, 0.12),       # L6, L7, L8
    (0.01, 0.12),                                   # c2
    (0.005, 0.040), (0.005, 0.023),                 # hsp, dsp
    (0.0, np.pi),                                   # theta_aux_fm
    (1.0, 8.0),                                     # gear_ratio
    (-np.pi, np.pi),                                # theta_offset
    (0.005, 0.025),                                 # Lpc
    (0.008, 0.030),                                 # Lpd
    (0.004, 0.030),                                 # Lac
    (0.0, np.deg2rad(170)),                         # BETA1 analitico
    (0.0, np.deg2rad(170)),                         # BETA2 analitico
    (0.0, np.pi)                                    # theta_aux_fd_fb
]

NOMBRES_PARAMETROS = [
    "Bancada1 (m)",     "Bancada2 (m)",
    "L1 (m)", "L2 (m)", "L3 (m)", "L4 (m)", "L5 (m)",
    "L6 (m)", "L7 (m)", "L8 (m)", "c2 (m)",
    "hsp (m)", "dsp (m)",
    "Theta Aux FM (rad)",
    "gear_ratio", "theta_offset (rad)",
    "Lpc (m)", "Lpd (m)", "Lac (m)",
    "BETA1 analitico (rad)", "BETA2 analitico (rad)",
    "Theta Aux FD fallback (rad)"
]


# ==============================================================================
# --- 7. ANGULOS GEOMETRICOS PARA CAD ---
# ==============================================================================
def reflect_about_line(point, line_start, line_end):
    """Refleja un punto sobre la recta que pasa por dos puntos."""
    d = line_end - line_start
    d = d / np.linalg.norm(d)
    v = point - line_start
    para = np.dot(v, d) * d
    perp = v - para
    return point - 2 * perp


def calcular_betas_geometricos(p, th_ref=0.0):
    """A partir del vector de parametros analiticos optimizado, calcula los
    angulos GEOMETRICOS BETA1_geom y BETA2_geom que se miden en el dibujo
    DORSAL real (los que entran al CAD)."""
    sim = run_kinematics(p, np.array([th_ref]))
    if sim is None:
        return None
    IFP = sim['ifp'][0]; IFD = sim['ifd'][0]; TIP = sim['tip'][0]

    (Bancada1, Bancada2, L1, L2, L3, L4, L5,
     L6, L7, L8, c2_par, hsp, dsp,
     theta_aux_fm, gear_ratio, theta_offset,
     Lpc, Lpd, Lac, BETA1, BETA2, _) = p

    r3 = Bancada1 / 2.0
    MCF = np.array([-r3, -Bancada2])

    # Recalcular CRK3 y ROK3 (analiticos, lado palmar)
    # Para esto necesitamos theta_fp, theta_fm, theta_fd en th_ref:
    th2 = th_ref
    th1 = (th2 / gear_ratio) + theta_offset
    theta14B = np.pi / 2.0
    c_arm = np.sqrt(hsp**2 + dsp**2)
    th4a = four_bar_open(L4, L5, c_arm, Bancada2, th1, theta14B)
    TH4a_deg = np.rad2deg(th4a)
    if TH4a_deg < 0:
        TH4a_deg += 360
    theta_fp = np.deg2rad(TH4a_deg + np.rad2deg(np.arctan2(hsp, dsp)))
    theta_fm = np.arctan2(IFD[1] - IFP[1], IFD[0] - IFP[0])
    theta_fd = np.arctan2(TIP[1] - IFD[1], TIP[0] - IFD[0])

    CRK3 = IFP + Lpc * np.array([np.cos(theta_fp + BETA1),
                                  np.sin(theta_fp + BETA1)])
    ROK3 = IFD + Lpd * np.array([np.cos(theta_fd + BETA2),
                                  np.sin(theta_fd + BETA2)])

    CRK3_dorsal = reflect_about_line(CRK3, IFP, IFD)
    ROK3_dorsal = reflect_about_line(ROK3, IFP, IFD)

    def ang(v_from, v_to):
        d = v_to - v_from
        return np.rad2deg(np.arctan2(d[1], d[0]))

    def ccw(a_from, a_to):
        return (a_to - a_from) % 360.0

    # BETA1_geom: de eje Fp (MCF -> IFP) a Lpc dorsal (IFP -> CRK3_dorsal)
    dir_fp = ang(MCF, IFP)
    dir_lpc_dorsal = ang(IFP, CRK3_dorsal)
    beta1_ccw = ccw(dir_fp, dir_lpc_dorsal)
    beta1_geom = beta1_ccw if beta1_ccw <= 180 else beta1_ccw - 360

    # BETA2_geom: de eje Fd (IFD -> TIP) a Lpd dorsal (IFD -> ROK3_dorsal)
    dir_fd = ang(IFD, TIP)
    dir_lpd_dorsal = ang(IFD, ROK3_dorsal)
    beta2_geom = ccw(dir_fd, dir_lpd_dorsal)

    return {
        'BETA1_geom_deg': beta1_geom,
        'BETA2_geom_deg': beta2_geom,
        'CRK3_dorsal': CRK3_dorsal,
        'ROK3_dorsal': ROK3_dorsal,
        'IFP': IFP, 'IFD': IFD, 'TIP': TIP,
        'theta_fp_deg': np.rad2deg(theta_fp),
        'theta_fm_deg': np.rad2deg(theta_fm),
        'theta_fd_deg': np.rad2deg(theta_fd),
    }


# ==============================================================================
# --- 8. OPTIMIZACION (Optuna + Differential Evolution) ---
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
    # --- VALIDACION: parametros equivalentes a CinematicaExoModificada.m ---
    p_matlab = [
        0.018, 0.020,                                  # Bancada1, Bancada2
        0.035, 0.049, 0.025, 0.020, 0.025,             # L1..L5
        0.055, 0.035, 0.052,                            # L6..L8
        0.04601,                                        # c2
        0.017, 0.018,                                   # hsp, dsp
        np.deg2rad(51.39),                             # theta_aux_fm
        2.0, np.deg2rad(109),                          # gear_ratio, theta_offset
        0.008, 0.018, 0.00886,                          # Lpc, Lpd, Lac
        np.deg2rad(40.0), np.deg2rad(110.0),            # BETA1, BETA2 analiticos
        np.deg2rad(38.78)                               # theta_aux_fd_fb
    ]
    print('\n>> VALIDACION con parametros equivalentes a CinematicaExoModificada.m')
    test_grid = np.linspace(0, np.deg2rad(132), 67)
    sim_test = run_kinematics(p_matlab, test_grid)
    if sim_test is None:
        print('   FALLO: cinematica no produce salida valida')
    else:
        # DIP relativo a partir de las posiciones
        ifp_t = sim_test['ifp']; ifd_t = sim_test['ifd']; tip_t = sim_test['tip']
        thfm_t = np.arctan2(ifd_t[:,1]-ifp_t[:,1], ifd_t[:,0]-ifp_t[:,0])
        thfd_t = np.arctan2(tip_t[:,1]-ifd_t[:,1], tip_t[:,0]-ifd_t[:,0])
        dip_t = np.rad2deg(np.unwrap(thfd_t - thfm_t))
        print(f'   DIP relativo: {dip_t[0]:.2f} -> {dip_t[-1]:.2f} deg '
              f'(excursion {dip_t.max()-dip_t.min():.2f})')
        diffs_t = np.diff(dip_t)
        mono = bool(np.all(diffs_t >= -1e-6) or np.all(diffs_t <= 1e-6))
        print(f'   Monotono: {"SI" if mono else "NO"}')
        info = calcular_betas_geometricos(p_matlab, th_ref=0.0)
        if info:
            print(f'   BETA1 analitico=40.000 -> BETA1_geom = {info["BETA1_geom_deg"]:8.3f} deg')
            print(f'   BETA2 analitico=110.000 -> BETA2_geom = {info["BETA2_geom_deg"]:8.3f} deg')

    # --- ETAPA 1: Optuna ---
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

    # --- ETAPA 2: Differential Evolution ---
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
    p_opt = resultado.x
    best_sim = run_kinematics(p_opt, theta_input)
    if best_sim is None:
        print('\n>> ADVERTENCIA: el resultado final no produce cinematica valida.')
        exit()

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

    print('\n====================================================')
    print('>> RESULTADOS')
    print('====================================================')
    print(f"   Error Global : {error_global_mm:.3f} mm")
    print(f"   Error IFP    : {err_ifp_mm:.3f} mm")
    print(f"   Error IFD    : {err_ifd_mm:.3f} mm")
    print(f"   Error Punta  : {err_tip_mm:.3f} mm\n")

    print('>> PARAMETROS DIMENSIONALES OPTIMIZADOS')
    for nombre, valor in zip(NOMBRES_PARAMETROS, p_opt):
        print(f"   {nombre:30s}: {valor:.6f}")

    print('\n>> PARAMETROS DE MONTAJE (transformacion rigida)')
    print(f'   Traslacion X : {t_opt[0]*1000:.2f} mm')
    print(f'   Traslacion Y : {t_opt[1]*1000:.2f} mm')
    print(f'   Rotacion Base: {angulo_montaje:.2f} deg')

    info = calcular_betas_geometricos(p_opt, th_ref=0.0)
    print('\n>> ANGULOS GEOMETRICOS PARA CAD (4B#3)')
    print('   (los que se miden en el dibujo DORSAL y entran al modelo CAD)')
    if info is not None:
        print(f"   BETA1 analitico = {np.rad2deg(p_opt[19]):8.3f} deg "
              f"-> BETA1_geom = {info['BETA1_geom_deg']:8.3f} deg  (eje Fp -> Lpc)")
        print(f"   BETA2 analitico = {np.rad2deg(p_opt[20]):8.3f} deg "
              f"-> BETA2_geom = {info['BETA2_geom_deg']:8.3f} deg  (eje Fd -> Lpd)")
        print(f"   Lpc = {p_opt[16]*1000:.3f} mm   "
              f"Lpd = {p_opt[17]*1000:.3f} mm   "
              f"Lac = {p_opt[18]*1000:.3f} mm")
    else:
        print("   (no se pudieron calcular: el 4B#3 no ensambla en theta=0)")
    print('====================================================\n')

    # ---- Grafica ----
    fig, ax = plt.subplots(figsize=(11, 8))
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
    ax.set_aspect('equal'); ax.grid(True, linestyle=':', alpha=0.7)
    ax.legend(loc='upper right', fontsize=9)
    ax.set_title(f'Sintesis con 4B#3 - Error Global: {error_global_mm:.3f} mm',
                 fontsize=14, fontweight='bold')
    ax.set_xlabel('Eje X (mm)', fontsize=12)
    ax.set_ylabel('Eje Y (mm)', fontsize=12)
    plt.tight_layout()
    plt.savefig('Resultados_Biofidelidad_4B3.png', dpi=150)
    plt.close()

    # ---- Guardado ----
    np.savetxt("Parametros_Optimizados_Mecanismo_4B3.txt", p_opt,
               header="Parametros optimizados (22 entradas, ver NOMBRES_PARAMETROS)")
    if info is not None:
        with open('Parametros_CAD_4B3.txt', 'w') as f:
            f.write("# Angulos geometricos del 4B#3 para usar en CAD (lado DORSAL)\n")
            f.write("# Estos son los angulos fisicos medidos con transportador en el dibujo,\n")
            f.write("# distintos a los analiticos (ver bloque de comentarios en exo2_1.py).\n\n")
            f.write(f"BETA1_geom_deg = {info['BETA1_geom_deg']:.4f}\n")
            f.write(f"BETA2_geom_deg = {info['BETA2_geom_deg']:.4f}\n")
            f.write(f"Lpc_mm = {p_opt[16]*1000:.4f}\n")
            f.write(f"Lpd_mm = {p_opt[17]*1000:.4f}\n")
            f.write(f"Lac_mm = {p_opt[18]*1000:.4f}\n")
            f.write(f"BETA1_analitico_deg = {np.rad2deg(p_opt[19]):.4f}\n")
            f.write(f"BETA2_analitico_deg = {np.rad2deg(p_opt[20]):.4f}\n")
        print(">> Archivo guardado: Parametros_CAD_4B3.txt (angulos para CAD)")

    print(">> Listo. Salidas:")
    print("   - Parametros_Optimizados_Mecanismo_4B3.txt (22 entradas)")
    print("   - Parametros_CAD_4B3.txt (angulos geometricos para CAD)")
    print("   - Resultados_Biofidelidad_4B3.png")
