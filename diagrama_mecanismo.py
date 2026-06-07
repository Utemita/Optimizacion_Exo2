"""
Diagrama de eslabones tipo INGENIERIA del exoesqueleto de dedo indice.

Reproduce las posiciones EXACTAS de cada pivote a partir de las ecuaciones de
CinematicaExoModificada.m en la posicion de referencia (THETA2=0).

Estilo de dibujo:
  - Lineas negras / gris oscuro para todos los eslabones.
  - Espesor grueso = cuerpos rigidos (falanges, marco).
  - Espesor medio = eslabones de los mecanismos.
  - Lineas a trazos = distancias virtuales / bancadas flotantes (c2, d2).
  - Circulos blancos vacios = articulaciones de revolucion.
  - Triangulos achurados = apoyos fijos (ground).
  - Arcos con flechas = angulos.
  - Etiquetas con nombre de parametro y valor numerico.
  - Mecanismos distribuidos sobre el dedo (dorsal=arriba).
  - Las falanges forman la linea inferior (palmar=abajo).

Un modelo CAD construido con las longitudes y angulos aqui indicados
reproduce EXACTAMENTE las trayectorias de CinematicaExoModificada.m.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import (Circle, Polygon, FancyArrowPatch, Arc,
                                FancyArrow, PathPatch)
from matplotlib.lines import Line2D
from matplotlib.path import Path

# ===================== PARAMETROS (dedo indice, mm/grados) =====================
Bancada1 = 18.0; Bancada2 = 20.0
Link1 = 35.0; Link2 = 49.0; Link3 = 25.0; Link4 = 20.0; Link5 = 25.0
Link6 = 55.0; Link7 = 35.0; Link8 = 52.0
c2 = 46.01
TETHA1inicial = 109.0; THETA14B = 90.0
hsp = 17.0; dsp = 18.0
fp = 49.0; fm = 26.0; fd = 24.0
THETAauxfm = 51.39; THETAauxfd = 38.78
reng = 2.0
# Tercer mecanismo 4B (DIP)
Lpc = 8.0; Lpd = 18.0; Lac = 8.86; BETA1 = 40.0; BETA2 = 110.0

# Variables derivadas
r1 = Link4; r2 = Link3; r3 = Bancada1 / 2.0; r4 = Link1; r5 = Link2
a = Link4; b = Link5; c = np.sqrt(hsp**2 + dsp**2); d = Bancada2
theta14B = np.deg2rad(THETA14B)
r1m2 = fp - 2*dsp; r2m2 = Link7; r3m2 = Link5 / 2.0; r4m2 = Link3; r5m2 = Link6
a2 = Link7; b2 = Link8; d2 = np.sqrt(hsp**2 + dsp**2)


# ===================== CINEMATICA =====================
def four_bar_open(a_, b_, c_, d_, th2_, th1_):
    """Resuelve 4-barras, configuracion abierta."""
    k1 = d_*np.cos(th1_) + a_*np.cos(th2_)
    k2 = d_*np.sin(th1_) + a_*np.sin(th2_)
    k3 = k1**2 + k2**2 + c_**2 - b_**2
    A1 = -k3 - 2*k1*c_
    B1 = 4*k2*c_
    C1 = 2*k1*c_ - k3
    disc = B1**2 - 4*A1*C1
    return 2*np.arctan((-B1 - np.sqrt(disc)) / (2*A1))


def _reflect_about_line(point, line_start, line_end):
    """Refleja un punto sobre una linea definida por dos puntos (preserva distancias)."""
    d_vec = line_end - line_start
    d_vec = d_vec / np.linalg.norm(d_vec)
    v = point - line_start
    para = np.dot(v, d_vec) * d_vec
    perp = v - para
    return point - 2 * perp


def compute_geometry(THETA2):
    """Calcula todas las posiciones de los puntos del mecanismo para un angulo dado."""
    th2 = np.deg2rad(THETA2)
    th1 = np.deg2rad(THETA2 / reng + TETHA1inicial)
    pts = {}

    # Marco fijo
    A = np.array([-r3, 0.0])
    B = np.array([r3, 0.0])
    MCF = np.array([-r3, -d])
    pts['A'] = A; pts['B'] = B; pts['MCF'] = MCF

    # Primer mecanismo de 5 barras
    M4 = A + r1*np.array([np.cos(th1), np.sin(th1)])
    J2 = B + r4*np.array([np.cos(th2), np.sin(th2)])
    e_ = (r1*np.sin(th1) - r4*np.sin(th2)) / (r4*np.cos(th2) - r1*np.cos(th1) + 2*r3)
    f_ = (2*(r1*r3*np.cos(th1) + r3*r4*np.cos(th2)) - r1**2 + r2**2 + r4**2 - r5**2) / \
         (2*(r4*np.cos(th2) - r1*np.cos(th1) + 2*r3))
    daux = e_**2 + 1
    g_ = 2*(e_*f_ - e_*r1*np.cos(th1) + e_*r3 - r1*np.sin(th1))
    h_ = f_**2 - 2*f_*(r1*np.cos(th1) - r3) - 2*r1*r3*np.cos(th1) + r1**2 + r3**2 - r2**2
    pyP = (-g_ + np.sqrt(g_**2 - 4*daux*h_)) / (2*daux)
    pxP = e_*pyP + f_
    Pp = np.array([pxP, pyP])
    pts['M4'] = M4; pts['J2'] = J2; pts['P'] = Pp

    # Primer mecanismo de 4 barras
    th4a = four_bar_open(a, b, c, d, th1, theta14B)
    TH4a = np.rad2deg(th4a)
    if TH4a < 0: TH4a += 360
    th4a = np.deg2rad(TH4a)
    thfp = np.deg2rad(TH4a + np.rad2deg(np.arctan2(hsp, dsp)))
    IFP = MCF + fp*np.array([np.cos(thfp), np.sin(thfp)])
    S1 = MCF + c*np.array([np.cos(th4a), np.sin(th4a)])
    thps2 = thfp - np.arctan2(hsp, fp - dsp)
    rs2 = np.sqrt(hsp**2 + (fp - dsp)**2)
    S2 = MCF + rs2*np.array([np.cos(thps2), np.sin(thps2)])
    pts['IFP'] = IFP; pts['S1'] = S1; pts['S2'] = S2; pts['thfp'] = thfp

    # Segundo mecanismo de 5 barras
    throll = np.arctan2(M4[1] - S1[1], M4[0] - S1[0])

    def ang_local(pt_to, pt_from):
        ang = np.rad2deg(np.arctan2(pt_to[1] - pt_from[1], pt_to[0] - pt_from[0]))
        if ang < 0: ang += 360
        return np.deg2rad(ang) - throll

    th1m2 = ang_local(S2, S1)
    th2m2 = ang_local(Pp, M4)
    em2 = (r1m2*np.sin(th1m2) - r4m2*np.sin(th2m2)) / \
          (r4m2*np.cos(th2m2) - r1m2*np.cos(th1m2) + 2*r3m2)
    fm2 = (2*(r1m2*r3m2*np.cos(th1m2) + r3m2*r4m2*np.cos(th2m2)) -
            r1m2**2 + r2m2**2 + r4m2**2 - r5m2**2) / \
           (2*(r4m2*np.cos(th2m2) - r1m2*np.cos(th1m2) + 2*r3m2))
    dauxm2 = em2**2 + 1
    gm2 = 2*(em2*fm2 - em2*r1m2*np.cos(th1m2) + em2*r3m2 - r1m2*np.sin(th1m2))
    hm2 = fm2**2 - 2*fm2*(r1m2*np.cos(th1m2) - r3m2) - \
           2*r1m2*r3m2*np.cos(th1m2) + r1m2**2 + r3m2**2 - r2m2**2
    pyP2 = (-gm2 + np.sqrt(gm2**2 - 4*dauxm2*hm2)) / (2*dauxm2)
    pxP2 = em2*pyP2 + fm2
    p2_ = np.hypot(pxP2, pyP2)
    th2p2 = np.arctan2(pyP2, pxP2)
    AUX = (S1 + M4) / 2.0
    P2 = p2_*np.array([np.cos(th2p2 + throll), np.sin(th2p2 + throll)]) + AUX
    pts['P2'] = P2

    # Segundo mecanismo de 4 barras
    th14B2 = np.arctan2(S2[1] - IFP[1], S2[0] - IFP[0])
    a24 = np.rad2deg(np.arctan2(P2[1] - S2[1], P2[0] - S2[0]))
    if a24 < 0: a24 += 360
    th24B2 = np.deg2rad(a24)
    th4am2 = four_bar_open(a2, b2, c2, d2, th24B2, th14B2)
    TH4am2 = np.rad2deg(th4am2)
    if TH4am2 < 0: TH4am2 += 360
    P3 = IFP + c2*np.array([np.cos(np.deg2rad(TH4am2)), np.sin(np.deg2rad(TH4am2))])
    thfm = np.deg2rad(TH4am2 + THETAauxfm)
    IFD = IFP + fm*np.array([np.cos(thfm), np.sin(thfm)])
    pts['P3'] = P3; pts['IFD'] = IFD; pts['thfm'] = thfm

    # Tercer mecanismo de 4 barras (DIP)
    thfp_ = pts['thfp']
    alpha1 = thfp_ + np.deg2rad(BETA1) - thfm
    CRK3 = IFP + Lpc*np.array([np.cos(thfp_ + np.deg2rad(BETA1)),
                                np.sin(thfp_ + np.deg2rad(BETA1))])
    Ax = Lpc*np.cos(alpha1); Ay = Lpc*np.sin(alpha1)
    Px = Ax - fm; Py = Ay; R = np.hypot(Px, Py)
    K = (Px**2 + Py**2 + Lpd**2 - Lac**2) / (2*Lpd)
    phi = np.arctan2(Py, Px)
    alpha2 = phi - np.arccos(np.clip(K / R, -1, 1))
    # Continuidad de angulo (unwrap respecto a llamada anterior)
    if hasattr(compute_geometry, '_alpha2_prev') and compute_geometry._alpha2_prev is not None:
        while alpha2 - compute_geometry._alpha2_prev > np.pi:
            alpha2 -= 2*np.pi
        while alpha2 - compute_geometry._alpha2_prev < -np.pi:
            alpha2 += 2*np.pi
    compute_geometry._alpha2_prev = alpha2
    thfd = thfm + (alpha2 - np.deg2rad(BETA2))
    ROK3 = IFD + Lpd*np.array([np.cos(thfd + np.deg2rad(BETA2)),
                                np.sin(thfd + np.deg2rad(BETA2))])
    TIP = IFD + fd*np.array([np.cos(thfd), np.sin(thfd)])

    # Reflejar CRK3 y ROK3 al lado DORSAL (sobre la linea IFP-IFD)
    # La solucion analitica los coloca del lado palmar; el montaje fisico
    # es dorsal. La reflexion preserva todas las longitudes.
    CRK3 = _reflect_about_line(CRK3, IFP, IFD)
    ROK3 = _reflect_about_line(ROK3, IFP, IFD)

    pts['CRK3'] = CRK3; pts['ROK3'] = ROK3; pts['TIP'] = TIP; pts['thfd'] = thfd

    return pts


# ===================== FUNCIONES DE DIBUJO =====================
def draw_link(ax, p1, p2, lw=2.0, ls='-', color='k', zorder=3):
    """Dibuja un eslabon como linea."""
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], ls=ls, color=color,
            lw=lw, zorder=zorder, solid_capstyle='round')


def draw_joint(ax, pos, radius=1.8, color='k', zorder=8):
    """Dibuja una articulacion de revolucion (circulo blanco con borde)."""
    circ = Circle(pos, radius, facecolor='white', edgecolor=color,
                  lw=1.5, zorder=zorder)
    ax.add_patch(circ)


def draw_ground(ax, pos, angle_deg=0, size=6.0, color='k', zorder=7):
    """Dibuja un soporte fijo (triangulo con achurado)."""
    ang = np.deg2rad(angle_deg)
    # Triangulo apuntando hacia abajo desde el punto
    base_dir = np.array([np.cos(ang), np.sin(ang)])
    perp = np.array([-np.sin(ang), np.cos(ang)])
    v0 = pos
    v1 = pos - size*base_dir + size*0.5*perp
    v2 = pos - size*base_dir - size*0.5*perp
    tri = Polygon([v0, v1, v2], closed=True, facecolor='none',
                  edgecolor=color, lw=1.5, zorder=zorder)
    ax.add_patch(tri)
    # Lineas de achurado
    base_center = pos - size*base_dir
    for i in range(4):
        t = (i + 0.5) / 4.0
        pt_on_base = v2 + t*(v1 - v2)
        hatch_end = pt_on_base - 2.5*base_dir
        ax.plot([pt_on_base[0], hatch_end[0]], [pt_on_base[1], hatch_end[1]],
                '-', color=color, lw=0.8, zorder=zorder - 1)


def draw_angle_arc(ax, center, angle_start_deg, angle_end_deg, radius=8.0,
                   label='', color='k', fontsize=7.5, label_offset=1.2):
    """Dibuja un arco con flecha para indicar un angulo."""
    arc = Arc(center, 2*radius, 2*radius, angle=0,
              theta1=angle_start_deg, theta2=angle_end_deg,
              color=color, lw=1.0, zorder=5)
    ax.add_patch(arc)
    # Flecha al final del arco
    mid_ang = np.deg2rad((angle_start_deg + angle_end_deg) / 2.0)
    label_pos = center + (radius + label_offset*3)*np.array([np.cos(mid_ang), np.sin(mid_ang)])
    if label:
        ax.annotate(label, label_pos, fontsize=fontsize, ha='center', va='center',
                    color=color, zorder=6)


def label_link(ax, p1, p2, text, offset=(0, 2.5), fontsize=7.5, color='k', ha='center'):
    """Coloca una etiqueta en el punto medio de un eslabon."""
    mid = (p1 + p2) / 2.0 + np.array(offset)
    ax.text(mid[0], mid[1], text, fontsize=fontsize, ha=ha, va='center',
            color=color, zorder=10,
            bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.85))


def label_point(ax, pos, text, dx=2.0, dy=2.0, fontsize=7.5, color='k',
                ha='left', va='bottom', bold=False):
    """Coloca una etiqueta junto a un punto."""
    weight = 'bold' if bold else 'normal'
    ax.annotate(text, (pos[0] + dx, pos[1] + dy), fontsize=fontsize,
                color=color, ha=ha, va=va, fontweight=weight, zorder=10)


# ===================== GENERAR DIAGRAMA =====================
P = compute_geometry(0.0)

fig, ax = plt.subplots(figsize=(18, 11))
ax.set_facecolor('white')

# Colores: todo en escala de grises/negro para estilo de ingenieria
BLK = '#000000'
DRK = '#333333'
GRY = '#555555'
LGR = '#777777'

# ---- FALANGES (cuerpos rigidos, lineas mas gruesas) ----
draw_link(ax, P['MCF'], P['IFP'], lw=6, color=DRK)
draw_link(ax, P['IFP'], P['IFD'], lw=6, color=DRK)
draw_link(ax, P['IFD'], P['TIP'], lw=6, color=DRK)

label_link(ax, P['MCF'], P['IFP'], 'Fp = %g mm' % fp, offset=(0, -5), color=DRK, fontsize=8)
label_link(ax, P['IFP'], P['IFD'], 'Fm = %g mm' % fm, offset=(0, -5), color=DRK, fontsize=8)
label_link(ax, P['IFD'], P['TIP'], 'Fd = %g mm' % fd, offset=(0, -5), color=DRK, fontsize=8)

# ---- MARCO FIJO / BANCADA ----
draw_link(ax, P['A'], P['B'], lw=3.5, color=BLK)
draw_link(ax, P['A'], P['MCF'], lw=3.5, color=BLK)

label_link(ax, P['A'], P['B'], 'B1 = 2r3 = %g mm' % Bancada1, offset=(0, 3), color=BLK, fontsize=8)
label_link(ax, P['A'], P['MCF'], 'B2 = d = %g mm' % Bancada2, offset=(-6, 0), color=BLK, fontsize=8)

# Soportes fijos (ground)
draw_ground(ax, P['A'], angle_deg=90, size=5.5)
draw_ground(ax, P['B'], angle_deg=90, size=5.5)
draw_ground(ax, P['MCF'], angle_deg=180, size=5.5)

# Engranes (indicacion esquematica)
ax.add_patch(Circle(P['B'], 5.0, fill=False, ec=BLK, lw=1.2, ls='--', zorder=2))
ax.add_patch(Circle(P['A'], 2.5, fill=False, ec=BLK, lw=1.2, ls='--', zorder=2))
ax.annotate('Motor\n(reng=%g)' % reng, (P['B'][0] + 7, P['B'][1] + 7),
            fontsize=8, color=BLK, ha='left', va='bottom',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=BLK, lw=0.8))

# ---- MECANISMO 5 BARRAS #1 ----
draw_link(ax, P['A'], P['M4'], lw=2.0, color=BLK)
draw_link(ax, P['M4'], P['P'], lw=2.0, color=BLK)
draw_link(ax, P['B'], P['J2'], lw=2.0, color=BLK)
draw_link(ax, P['J2'], P['P'], lw=2.0, color=BLK)

label_link(ax, P['A'], P['M4'], 'L4 = %g' % Link4, offset=(-3, 2), fontsize=7.5)
label_link(ax, P['M4'], P['P'], 'L3 = %g' % Link3, offset=(0, 2.5), fontsize=7.5)
label_link(ax, P['B'], P['J2'], 'L1 = %g' % Link1, offset=(3, 1), fontsize=7.5)
label_link(ax, P['J2'], P['P'], 'L2 = %g' % Link2, offset=(2, 2), fontsize=7.5)

# ---- MECANISMO 4 BARRAS #1 (soporte falange proximal) ----
draw_link(ax, P['M4'], P['S1'], lw=2.0, color=BLK)
draw_link(ax, P['MCF'], P['S1'], lw=2.0, color=BLK)

label_link(ax, P['M4'], P['S1'], 'L5 = %g' % Link5, offset=(2, 2), fontsize=7.5)
cval = np.sqrt(hsp**2 + dsp**2)
label_link(ax, P['MCF'], P['S1'], 'c = %.1f' % cval, offset=(-3, -2), fontsize=7.5)

# Soportes S1, S2 sobre la falange proximal
# S1 y S2 son pines de articulacion sobre la falange con altura hsp y distancia dsp
draw_link(ax, P['S1'], P['S2'], lw=1.5, ls='--', color=GRY)
label_link(ax, P['S1'], P['S2'], 'fp-2dsp = %g' % (fp - 2*dsp), offset=(0, 3), fontsize=7, color=GRY)

# ---- MECANISMO 5 BARRAS #2 ----
draw_link(ax, P['S2'], P['P2'], lw=2.0, color=BLK)
draw_link(ax, P['P'], P['P2'], lw=2.0, color=BLK)

label_link(ax, P['S2'], P['P2'], 'L7 = %g' % Link7, offset=(2, 2.5), fontsize=7.5)
label_link(ax, P['P'], P['P2'], 'L6 = %g' % Link6, offset=(0, 3), fontsize=7.5)

# ---- MECANISMO 4 BARRAS #2 (driver medial) ----
draw_link(ax, P['P2'], P['P3'], lw=2.0, color=BLK)
draw_link(ax, P['IFP'], P['P3'], lw=2.0, ls='--', color=BLK)
draw_link(ax, P['IFP'], P['S2'], lw=1.5, ls='--', color=GRY)

label_link(ax, P['P2'], P['P3'], 'L8 = %g' % Link8, offset=(2, 2.5), fontsize=7.5)
label_link(ax, P['IFP'], P['P3'], 'c2 = %.2f' % c2, offset=(2, 2), fontsize=7.5)
label_link(ax, P['IFP'], P['S2'], 'd2 = %.1f (bancada)' % d2, offset=(0, -2.5), fontsize=7, color=GRY)

# ---- MECANISMO 4 BARRAS #3 (DIP) ----
draw_link(ax, P['IFP'], P['CRK3'], lw=2.5, color=BLK)
draw_link(ax, P['CRK3'], P['ROK3'], lw=2.5, color=BLK)
draw_link(ax, P['IFD'], P['ROK3'], lw=2.5, color=BLK)
# La bancada del 4B#3 es la falange medial (IFP -> IFD), ya dibujada

label_link(ax, P['IFP'], P['CRK3'], 'Lpc = %g' % Lpc, offset=(-1, 2.5), fontsize=7.5)
label_link(ax, P['CRK3'], P['ROK3'], 'Lac = %.2f' % Lac, offset=(0, 2.5), fontsize=7.5)
label_link(ax, P['IFD'], P['ROK3'], 'Lpd = %g' % Lpd, offset=(1, 2.5), fontsize=7.5)

# ---- ARTICULACIONES (circulos blancos) ----
joints_main = {
    'A': ('A', -3, 4),
    'B': ('B', 2, 4),
    'MCF': ('MCF', -3, -4),
    'M4': ('M4', 2, 2.5),
    'J2': ('J2', 2, -3),
    'P': ('P', 2, 2.5),
    'S1': ('S1', -3, 3),
    'S2': ('S2', 2, 3),
    'P2': ('P2', 2, 2.5),
    'P3': ('P3', 2, 2.5),
    'IFP': ('IFP', 2, -4),
    'IFD': ('IFD', 2, -4),
    'CRK3': ('CRK3', -2, 2.5),
    'ROK3': ('ROK3', 2, 2.5),
    'TIP': ('TIP', 2, -3),
}
for key, (lbl, dx, dy) in joints_main.items():
    draw_joint(ax, P[key], radius=1.6)
    label_point(ax, P[key], lbl, dx=dx, dy=dy, fontsize=7, bold=True)

# ---- ANGULOS ----
# BETA1: angulo de montaje de la manivela respecto a la proximal en IFP
thfp_deg = np.rad2deg(P['thfp'])
draw_angle_arc(ax, P['IFP'], thfp_deg, thfp_deg + BETA1, radius=10, color=GRY)
label_pt_beta1 = P['IFP'] + 14*np.array([
    np.cos(np.deg2rad(thfp_deg + BETA1/2)),
    np.sin(np.deg2rad(thfp_deg + BETA1/2))])
ax.text(label_pt_beta1[0], label_pt_beta1[1], 'BETA1\n=%g deg' % BETA1,
        fontsize=6.5, ha='center', va='center', color=GRY)

# BETA2: angulo de montaje del balancin respecto a la distal en IFD
thfd_deg = np.rad2deg(P['thfd'])
draw_angle_arc(ax, P['IFD'], thfd_deg, thfd_deg + BETA2, radius=10, color=GRY)
label_pt_beta2 = P['IFD'] + 14*np.array([
    np.cos(np.deg2rad(thfd_deg + BETA2/2)),
    np.sin(np.deg2rad(thfd_deg + BETA2/2))])
ax.text(label_pt_beta2[0], label_pt_beta2[1], 'BETA2\n=%g deg' % BETA2,
        fontsize=6.5, ha='center', va='center', color=GRY)

# theta_1_inicial
draw_angle_arc(ax, P['A'], 0, TETHA1inicial, radius=7, color=LGR)
lbl_th1 = P['A'] + 10*np.array([np.cos(np.deg2rad(TETHA1inicial/2)),
                                  np.sin(np.deg2rad(TETHA1inicial/2))])
ax.text(lbl_th1[0], lbl_th1[1], 'theta1_ini\n=%g deg' % TETHA1inicial,
        fontsize=6, ha='center', va='center', color=LGR)

# THETAauxfm
th4am2_deg = np.rad2deg(np.arctan2(P['P3'][1] - P['IFP'][1], P['P3'][0] - P['IFP'][0]))
thfm_deg = np.rad2deg(P['thfm'])
draw_angle_arc(ax, P['IFP'], th4am2_deg, thfm_deg, radius=12, color=LGR)
lbl_auxfm = P['IFP'] + 16*np.array([np.cos(np.deg2rad((th4am2_deg + thfm_deg)/2)),
                                      np.sin(np.deg2rad((th4am2_deg + thfm_deg)/2))])
ax.text(lbl_auxfm[0], lbl_auxfm[1], 'THETAauxfm\n=%.2f deg' % THETAauxfm,
        fontsize=6, ha='center', va='center', color=LGR)

# ---- ANOTACIONES DE SOPORTE ----
# hsp y dsp
ax.annotate('hsp = %g mm\ndsp = %g mm' % (hsp, dsp),
            xy=P['S1'], xytext=(P['S1'][0] - 12, P['S1'][1] + 8),
            fontsize=7, color=GRY, ha='center',
            arrowprops=dict(arrowstyle='->', color=GRY, lw=0.8),
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=GRY, lw=0.6))

# ---- LEYENDA (estilo de ingenieria) ----
legend_elements = [
    Line2D([0], [0], color=DRK, lw=6, label='Cuerpos rigidos (falanges)'),
    Line2D([0], [0], color=BLK, lw=3.5, label='Marco fijo (bancada)'),
    Line2D([0], [0], color=BLK, lw=2.0, label='Eslabones de mecanismos'),
    Line2D([0], [0], color=BLK, lw=2.0, ls='--', label='Distancias virtuales / bancadas flotantes'),
    Line2D([0], [0], color='white', marker='o', markeredgecolor=BLK,
           markerfacecolor='white', markersize=8, lw=0, label='Articulacion de revolucion'),
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=9,
          framealpha=0.95, edgecolor=BLK)

# ---- CUADRO DE INFORMACION del 4B#3 ----
info_text = (
    '4-Barras #3 (DIP) - DORSAL:\n'
    '  Posicion: lado DORSAL (reflejado del analitico)\n'
    '  Bancada = Fm = %g mm (IFP a IFD)\n'
    '  Manivela (Lpc=%g) rigida a Fp, pivote en IFP\n'
    '  Acoplador (Lac=%.2f) CRK3-ROK3\n'
    '  Balancin (Lpd=%g) rigida a Fd, pivote en IFD\n'
    '  BETA1=%g deg, BETA2=%g deg\n'
    '  Entrada: rotacion relativa Fp vs Fm (PIP)' %
    (fm, Lpc, Lac, Lpd, BETA1, BETA2)
)
ax.text(0.99, 0.02, info_text, transform=ax.transAxes,
        fontsize=8, va='bottom', ha='right', color=BLK,
        bbox=dict(boxstyle='round,pad=0.4', fc='#f8f8f8', ec=BLK, lw=1.0),
        fontfamily='monospace', zorder=10)

# ---- CUADRO DE PARAMETROS COMPLETOS ----
param_text = (
    'PARAMETROS COMPLETOS (mm, deg)\n'
    '------------------------------\n'
    'Bancada1 = %g   Bancada2 = %g\n'
    'Link1 = %g  Link2 = %g  Link3 = %g\n'
    'Link4 = %g  Link5 = %g  Link6 = %g\n'
    'Link7 = %g  Link8 = %g  c2 = %.2f\n'
    'hsp = %g  dsp = %g  fp = %g\n'
    'fm = %g  fd = %g\n'
    'THETA1_ini = %g  THETA14B = %g\n'
    'THETAauxfm = %.2f  reng = %g\n'
    'Lpc = %g  Lpd = %g  Lac = %.2f\n'
    'BETA1 = %g  BETA2 = %g' %
    (Bancada1, Bancada2, Link1, Link2, Link3, Link4, Link5, Link6,
     Link7, Link8, c2, hsp, dsp, fp, fm, fd,
     TETHA1inicial, THETA14B, THETAauxfm, reng,
     Lpc, Lpd, Lac, BETA1, BETA2)
)
ax.text(0.01, 0.98, param_text, transform=ax.transAxes,
        fontsize=7.5, va='top', ha='left', color=BLK,
        bbox=dict(boxstyle='round,pad=0.4', fc='#f8f8ff', ec=BLK, lw=1.0),
        fontfamily='monospace', zorder=10)

# ---- TITULO Y CONFIGURACION ----
ax.set_aspect('equal')
ax.grid(True, ls=':', alpha=0.4, color='#999999')
ax.set_title('Diagrama de Eslabones - Exoesqueleto de Dedo Indice\n'
             'Posicion de referencia: THETA2 = 0 deg (cotas en mm)\n'
             'Construir segun este diagrama reproduce las trayectorias de CinematicaExoModificada.m',
             fontsize=12, fontweight='bold', color=BLK)
ax.set_xlabel('X (mm)', fontsize=10)
ax.set_ylabel('Y (mm)', fontsize=10)

plt.tight_layout()
plt.savefig('diagrama_mecanismo_completo.png', dpi=180, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()

# ===================== VERIFICACION DE LONGITUDES =====================
def dist(a_, b_): return np.linalg.norm(P[a_] - P[b_])

checks = [
    ("|A-M4| = Link4",              dist('A', 'M4'),    Link4),
    ("|M4-P| = Link3",              dist('M4', 'P'),    Link3),
    ("|B-J2| = Link1",              dist('B', 'J2'),    Link1),
    ("|J2-P| = Link2",              dist('J2', 'P'),    Link2),
    ("|M4-S1| = Link5",             dist('M4', 'S1'),   Link5),
    ("|MCF-S1| = c",                dist('MCF', 'S1'),  np.sqrt(hsp**2 + dsp**2)),
    ("|S1-S2| = fp-2dsp",           dist('S1', 'S2'),   fp - 2*dsp),
    ("|S2-P2| = Link7",             dist('S2', 'P2'),   Link7),
    ("|P-P2| = Link6",              dist('P', 'P2'),    Link6),
    ("|P2-P3| = Link8",             dist('P2', 'P3'),   Link8),
    ("|IFP-P3| = c2",               dist('IFP', 'P3'),  c2),
    ("|MCF-IFP| = fp",              dist('MCF', 'IFP'), fp),
    ("|IFP-IFD| = fm",              dist('IFP', 'IFD'), fm),
    ("|IFD-TIP| = fd",              dist('IFD', 'TIP'), fd),
    ("|IFP-CRK3| = Lpc",            dist('IFP', 'CRK3'), Lpc),
    ("|CRK3-ROK3| = Lac",           dist('CRK3', 'ROK3'), Lac),
    ("|IFD-ROK3| = Lpd",            dist('IFD', 'ROK3'), Lpd),
]
print("=== VERIFICACION DE LONGITUDES (mm) ===")
all_ok = True
for name, val, exp in checks:
    ok = abs(val - exp) < 1e-4
    all_ok = all_ok and ok
    print(f"  [{'OK' if ok else 'XX'}] {name:30s}  calc={val:8.3f}  esperado={exp:8.3f}")
print(f"\n{'TODAS CORRECTAS' if all_ok else 'HAY DISCREPANCIAS - REVISAR'}")

# Verificacion DIP
print("\n=== VERIFICACION DIP (THETA2 = 0..132 deg) ===")
compute_geometry._alpha2_prev = None  # Reset unwrap state
dips = []
for T in np.linspace(0, 132, 67):
    PT = compute_geometry(T)
    dips.append(np.rad2deg(PT['thfd'] - PT['thfm']))
dips = np.array(dips)
diffs = np.diff(dips)
mono = bool(np.all(diffs >= -1e-9) or np.all(diffs <= 1e-9))
print(f"  DIP relativo: {dips[0]:.2f} -> {dips[-1]:.2f} deg  (excursion {dips.max()-dips.min():.2f} deg)")
print(f"  Monotono: {'SI' if mono else 'NO'}")

# Verificacion DORSAL: CRK3 y ROK3 deben estar del lado dorsal de la falange medial
print("\n=== VERIFICACION DORSAL (CRK3, ROK3 arriba de la linea IFP-IFD) ===")
d_phal = P['IFD'] - P['IFP']
d_norm = d_phal / np.linalg.norm(d_phal)
dorsal_dir = np.array([d_norm[1], -d_norm[0]])  # CW rotation = dorsal (arriba)
dot_crk3 = np.dot(P['CRK3'] - P['IFP'], dorsal_dir)
dot_rok3 = np.dot(P['ROK3'] - P['IFP'], dorsal_dir)
print(f"  CRK3 dot dorsal = {dot_crk3:.3f}  ({'DORSAL OK' if dot_crk3 > 0 else 'PALMAR - ERROR'})")
print(f"  ROK3 dot dorsal = {dot_rok3:.3f}  ({'DORSAL OK' if dot_rok3 > 0 else 'PALMAR - ERROR'})")

print(f"\nDiagrama guardado: diagrama_mecanismo_completo.png")
