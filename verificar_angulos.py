"""Verifica los angulos geometricos REALES de Lpc y Lpd despues de la reflexion."""
import sys
import numpy as np
sys.path.insert(0, '/projects/sandbox/Optimizacion_Exo2')
from diagrama_mecanismo import compute_geometry, BETA1, BETA2

P = compute_geometry(0.0)

# Vectores
IFP = P['IFP']; IFD = P['IFD']; TIP = P['TIP']
MCF = P['MCF']
CRK3 = P['CRK3']; ROK3 = P['ROK3']
thfp = P['thfp']; thfm = P['thfm']; thfd = P['thfd']

print("=== POSICIONES (mm) en THETA2=0 ===")
print(f"  MCF  = {MCF}")
print(f"  IFP  = {IFP}")
print(f"  IFD  = {IFD}")
print(f"  TIP  = {TIP}")
print(f"  CRK3 = {CRK3}  (extremo de Lpc, despues de reflejar a dorsal)")
print(f"  ROK3 = {ROK3}  (extremo de Lpd, despues de reflejar a dorsal)")

print("\n=== ANGULOS DE FALANGES (deg, global) ===")
print(f"  thfp (Fp: MCF->IFP) = {np.rad2deg(thfp):.2f}")
print(f"  thfm (Fm: IFP->IFD) = {np.rad2deg(thfm):.2f}")
print(f"  thfd (Fd: IFD->TIP) = {np.rad2deg(thfd):.2f}")

# Direccion de Fp (de MCF hacia IFP) y direccion del poste Lpc real (de IFP hacia CRK3 dorsal)
def angle_deg(v_from, v_to):
    d = v_to - v_from
    return np.rad2deg(np.arctan2(d[1], d[0]))

dir_fp = angle_deg(MCF, IFP)
dir_fm = angle_deg(IFP, IFD)
dir_fd = angle_deg(IFD, TIP)
dir_lpc_real = angle_deg(IFP, CRK3)
dir_lpd_real = angle_deg(IFD, ROK3)

print("\n=== DIRECCIONES (deg, global) ===")
print(f"  Fp (MCF->IFP)  = {dir_fp:.2f}")
print(f"  Fm (IFP->IFD)  = {dir_fm:.2f}")
print(f"  Fd (IFD->TIP)  = {dir_fd:.2f}")
print(f"  Lpc dorsal (IFP->CRK3 reflejado) = {dir_lpc_real:.2f}")
print(f"  Lpd dorsal (IFD->ROK3 reflejado) = {dir_lpd_real:.2f}")

# Angulos GEOMETRICOS (CCW desde direccion de la falange hasta la direccion del poste)
def ang_ccw(a_from, a_to):
    """Angulo CCW desde a_from hasta a_to, en [0, 360)."""
    return (a_to - a_from) % 360

beta1_geom = ang_ccw(dir_fp, dir_lpc_real)
beta2_geom = ang_ccw(dir_fd, dir_lpd_real)

print("\n=== ANGULOS GEOMETRICOS (lo que se mide en el diagrama) ===")
print(f"  BETA1 geometrico (de Fp a Lpc dorsal, CCW) = {beta1_geom:.2f} deg")
print(f"  BETA2 geometrico (de Fd a Lpd dorsal, CCW) = {beta2_geom:.2f} deg")

print("\n=== COMPARACION ===")
print(f"  BETA1 en codigo = {BETA1} deg   |   BETA1 fisico/geometrico = {beta1_geom:.2f} deg")
print(f"  BETA2 en codigo = {BETA2} deg   |   BETA2 fisico/geometrico = {beta2_geom:.2f} deg")

print("\n=== CHEQUEO LADO DORSAL ===")
# La normal al eje del dedo (apuntando dorsal) es Fm rotado +90 grados (CCW visto desde +Z)
# Como el dedo va de MCF hacia abajo a la izquierda, dorsal es hacia arriba
n_fm = np.array([-np.sin(np.deg2rad(dir_fm)), np.cos(np.deg2rad(dir_fm))])
print(f"  Normal dorsal Fm = {n_fm}  (debe apuntar +Y)")
print(f"  CRK3 dot normal = {np.dot(CRK3 - IFP, n_fm):.3f}  (>0 dorsal)")
print(f"  ROK3 dot normal = {np.dot(ROK3 - IFP, n_fm):.3f}  (>0 dorsal)")
