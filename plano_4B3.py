"""
Plano tecnico (vista lateral cotada) del tercer mecanismo de 4 barras (4B#3).

Genera dos figuras:
  1. plano_4B3_general.png    - Vista lateral del 4B#3 montado sobre las falanges
                                 con BETA1_geom, BETA2_geom, Lpc, Lac, Lpd cotados.
  2. plano_4B3_piezas.png     - Vistas de las piezas individuales con todas las
                                 dimensiones para fabricacion (acoplador y postes).
"""
import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch, Arc, Polygon
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from diagrama_mecanismo import (
    compute_geometry, Lpc, Lpd, Lac, BETA1, BETA2, fm, fd, fp, hsp,
)

# ===================== FIGURA 1: VISTA GENERAL COTADA DEL 4B#3 =====================
P = compute_geometry(0.0)
IFP = P['IFP']; IFD = P['IFD']; TIP = P['TIP']
MCF = P['MCF']
CRK3 = P['CRK3']; ROK3 = P['ROK3']

# Calcular angulos geometricos reales
def _ang(v_from, v_to):
    d = v_to - v_from
    return np.rad2deg(np.arctan2(d[1], d[0]))


def _ccw(a_from, a_to):
    return (a_to - a_from) % 360


dir_fp = _ang(MCF, IFP)
dir_fm = _ang(IFP, IFD)
dir_fd = _ang(IFD, TIP)
dir_lpc = _ang(IFP, CRK3)
dir_lpd = _ang(IFD, ROK3)

beta1_geom_ccw = _ccw(dir_fp, dir_lpc)
beta1_geom = beta1_geom_ccw if beta1_geom_ccw <= 180 else beta1_geom_ccw - 360
beta2_geom = _ccw(dir_fd, dir_lpd)


def draw_dim_line(ax, p1, p2, label, offset_dir=None, offset_mag=4, fontsize=8,
                   color='#1565c0', label_offset_extra=0, label_above=True):
    """Dibuja una cota (linea con flechas dobles y etiqueta) entre p1 y p2."""
    p1 = np.array(p1); p2 = np.array(p2)
    if offset_dir is None:
        # Direccion perpendicular automatica
        v = p2 - p1
        n = np.array([-v[1], v[0]]) / np.linalg.norm(v)
    else:
        n = np.array(offset_dir) / np.linalg.norm(offset_dir)

    p1d = p1 + n * offset_mag
    p2d = p2 + n * offset_mag
    # Lineas de extension
    ax.plot([p1[0], p1d[0]], [p1[1], p1d[1]], color=color, lw=0.6)
    ax.plot([p2[0], p2d[0]], [p2[1], p2d[1]], color=color, lw=0.6)
    # Linea de cota con flechas
    arrow = FancyArrowPatch(p1d, p2d, arrowstyle='<|-|>', color=color, lw=1.0,
                             mutation_scale=8)
    ax.add_patch(arrow)
    # Etiqueta
    mid = (p1d + p2d) / 2 + n * label_offset_extra
    ax.text(mid[0], mid[1], label, color=color, fontsize=fontsize,
            ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.9))


def draw_angle_arc_dim(ax, center, ang_start_deg, ang_end_deg, radius, label, color='#c62828',
                       fontsize=8):
    arc = Arc(center, 2*radius, 2*radius, angle=0,
              theta1=ang_start_deg, theta2=ang_end_deg, color=color, lw=1.0)
    ax.add_patch(arc)
    mid = np.deg2rad((ang_start_deg + ang_end_deg) / 2)
    label_pos = np.array(center) + (radius + 4) * np.array([np.cos(mid), np.sin(mid)])
    ax.text(label_pos[0], label_pos[1], label, color=color, fontsize=fontsize,
            ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec=color, lw=0.6))


# Figura
fig1, ax1 = plt.subplots(figsize=(14, 10))
ax1.set_facecolor('white')

DRK = '#333333'
ORG = '#cc6600'

# Falanges (simplificadas, solo Fp final, Fm completa, Fd completa)
ax1.plot([MCF[0], IFP[0]], [MCF[1], IFP[1]], color=DRK, lw=5, alpha=0.4)
ax1.plot([IFP[0], IFD[0]], [IFP[1], IFD[1]], color=DRK, lw=5)
ax1.plot([IFD[0], TIP[0]], [IFD[1], TIP[1]], color=DRK, lw=5)

# Articulaciones
for pt, lbl, dx, dy in [(IFP, 'IFP', 2, -4), (IFD, 'IFD', 2, -4),
                         (CRK3, 'CRK3', -2, 3), (ROK3, 'ROK3', 2, 3),
                         (TIP, 'TIP', 2, -3)]:
    ax1.add_patch(Circle(pt, 1.2, fc='white', ec='black', lw=1.2, zorder=10))
    ax1.text(pt[0] + dx, pt[1] + dy, lbl, fontsize=9, fontweight='bold')

# 4B#3
ax1.plot([IFP[0], CRK3[0]], [IFP[1], CRK3[1]], color='black', lw=2.0)
ax1.plot([CRK3[0], ROK3[0]], [CRK3[1], ROK3[1]], color='black', lw=2.0)
ax1.plot([IFD[0], ROK3[0]], [IFD[1], ROK3[1]], color='black', lw=2.0)

# Cotas de longitudes
draw_dim_line(ax1, IFP, CRK3, 'Lpc = %g mm' % Lpc, offset_mag=5, label_above=True)
draw_dim_line(ax1, CRK3, ROK3, 'Lac = %.2f mm' % Lac, offset_mag=4)
draw_dim_line(ax1, IFD, ROK3, 'Lpd = %g mm' % Lpd, offset_mag=5)
draw_dim_line(ax1, IFP, IFD, 'Fm = %g mm (bancada del 4B#3)' % fm, offset_mag=-7,
               color='#5d4037', label_above=False)

# Angulos BETA1_geom y BETA2_geom
draw_angle_arc_dim(ax1, IFP, dir_lpc, dir_fp, radius=10,
                    label='BETA1_geom = %.1f deg' % beta1_geom)
draw_angle_arc_dim(ax1, IFD, dir_fd, dir_fd + beta2_geom, radius=12,
                    label='BETA2_geom = %.1f deg' % beta2_geom)

# Direcciones de referencia (Fp, Fd) como lineas a trazos
ax1.plot([IFP[0], IFP[0] + 18*np.cos(np.deg2rad(dir_fp))],
         [IFP[1], IFP[1] + 18*np.sin(np.deg2rad(dir_fp))],
         color='#999', ls=':', lw=0.8)
ax1.plot([IFD[0], IFD[0] + 18*np.cos(np.deg2rad(dir_fd))],
         [IFD[1], IFD[1] + 18*np.sin(np.deg2rad(dir_fd))],
         color='#999', ls=':', lw=0.8)
ax1.text(IFP[0] + 18*np.cos(np.deg2rad(dir_fp)),
         IFP[1] + 18*np.sin(np.deg2rad(dir_fp)) + 2,
         'eje Fp (hacia MCF)', fontsize=7, color='#666')
ax1.text(IFD[0] + 18*np.cos(np.deg2rad(dir_fd)) + 1,
         IFD[1] + 18*np.sin(np.deg2rad(dir_fd)),
         'eje Fd (hacia TIP)', fontsize=7, color='#666')

ax1.set_aspect('equal')
ax1.grid(True, ls=':', alpha=0.4)
ax1.set_title('Plano tecnico - 4B#3 (DIP) montaje DORSAL\n'
              'Posicion de referencia THETA2 = 0 deg', fontsize=11, fontweight='bold')
ax1.set_xlabel('X [mm]')
ax1.set_ylabel('Y [mm]')

# Cuadro de notas
notes = (
    'NOTAS DE CONSTRUCCION:\n'
    '----------------------------------------\n'
    '1. Postes Lpc y Lpd son OREJAS RIGIDAS\n'
    '   atornilladas a la falange respectiva.\n'
    '   No giran sobre la falange: la falange\n'
    '   y su poste se mueven como un solo cuerpo.\n'
    '2. CRK3 y ROK3 son articulaciones de\n'
    '   revolucion (pasadores).\n'
    '3. Acoplador Lac libre, conectado por\n'
    '   pasadores en CRK3 y ROK3.\n'
    '4. Mecanismo del lado DORSAL del dedo.\n'
    '   Ningun eslabon invade lado palmar.\n'
    '5. BETA1_geom y BETA2_geom son los angulos\n'
    '   medidos en este dibujo.'
)
ax1.text(0.01, 0.98, notes, transform=ax1.transAxes,
         fontsize=8, va='top', ha='left',
         bbox=dict(boxstyle='round,pad=0.4', fc='#fffde7', ec='#999', lw=0.8),
         fontfamily='monospace')

plt.tight_layout()
plt.savefig('plano_4B3_general.png', dpi=200, bbox_inches='tight',
            facecolor='white')
plt.close()

print("Figura 1 guardada: plano_4B3_general.png")


# ===================== FIGURA 2: PIEZAS INDIVIDUALES COTADAS =====================
fig2 = plt.figure(figsize=(14, 11))

# Constantes de fabricacion (las mismas que en el script de STEP)
T_LINK = 1.5      # espesor de eslabones/orejas
W_LINK = 4.0      # ancho de orejas
END_PAD = 1.5     # padding desde borde de agujero al final
D_PIN = 1.5       # diametro del pin
D_HOLE = D_PIN + 0.20   # diametro del agujero
BASE_LEN = 6.0
BASE_WID = 4.0
BASE_THK = 1.0
SCREW_DIA = 1.5

# ---- Subplot 1: Acoplador ----
ax2 = plt.subplot(3, 2, 1)
L_total = Lac + 2 * (D_HOLE / 2 + END_PAD)
# Cuerpo
rect = Rectangle((-L_total/2, -W_LINK/2), L_total, W_LINK,
                  fc='#bbdefb', ec='black', lw=1.5)
ax2.add_patch(rect)
# Agujeros
ax2.add_patch(Circle((-Lac/2, 0), D_HOLE/2, fc='white', ec='black', lw=1.2))
ax2.add_patch(Circle((Lac/2, 0), D_HOLE/2, fc='white', ec='black', lw=1.2))
# Cotas
draw_dim_line(ax2, (-Lac/2, 0), (Lac/2, 0), 'Lac = %.2f' % Lac, offset_mag=-3.5,
               fontsize=8)
draw_dim_line(ax2, (-L_total/2, -W_LINK/2), (L_total/2, -W_LINK/2),
               'L_tot = %.2f' % L_total, offset_mag=-3.5, fontsize=7,
               color='#666')
draw_dim_line(ax2, (-L_total/2 - 0.3, -W_LINK/2), (-L_total/2 - 0.3, W_LINK/2),
               'W=%g' % W_LINK, offset_mag=-3.5, fontsize=7, color='#666')
ax2.annotate('agujero %g' % D_HOLE,
             xy=(-Lac/2, 0), xytext=(-Lac/2 - 2, W_LINK + 1),
             fontsize=7, color='#c62828',
             arrowprops=dict(arrowstyle='->', color='#c62828', lw=0.6))

ax2.set_xlim(-L_total/2 - 5, L_total/2 + 5)
ax2.set_ylim(-W_LINK - 2, W_LINK + 2.5)
ax2.set_aspect('equal')
ax2.grid(True, ls=':', alpha=0.3)
ax2.set_title('ACOPLADOR Lac\n(vista en planta)', fontsize=10, fontweight='bold')
ax2.set_xlabel('mm'); ax2.set_ylabel('mm')

# ---- Subplot 1b: Acoplador vista lateral (espesor) ----
ax2b = plt.subplot(3, 2, 2)
rect_side = Rectangle((-L_total/2, 0), L_total, T_LINK,
                        fc='#bbdefb', ec='black', lw=1.5)
ax2b.add_patch(rect_side)
draw_dim_line(ax2b, (-L_total/2 - 0.3, 0), (-L_total/2 - 0.3, T_LINK),
               'T = %g' % T_LINK, offset_mag=-2, fontsize=7, color='#666')
ax2b.set_xlim(-L_total/2 - 4, L_total/2 + 4)
ax2b.set_ylim(-2, 6)
ax2b.set_aspect('equal')
ax2b.grid(True, ls=':', alpha=0.3)
ax2b.set_title('Acoplador Lac (vista lateral - espesor)', fontsize=10)
ax2b.set_xlabel('mm'); ax2b.set_ylabel('mm')


# ---- Subplot 2: Poste Lpc ----
ax3 = plt.subplot(3, 2, 3)
body_len = Lpc + D_HOLE / 2 + END_PAD
# Base de montaje (rectangulo en -X)
base_x_min = -BASE_LEN
base_y_min = -BASE_WID/2
ax3.add_patch(Rectangle((base_x_min, base_y_min), BASE_LEN, BASE_WID,
                          fc='#c8e6c9', ec='black', lw=1.5))
# Cuerpo del poste (rectangulo desde 0 hasta body_len)
ax3.add_patch(Rectangle((0, -W_LINK/2), body_len, W_LINK,
                          fc='#c8e6c9', ec='black', lw=1.5))
# Agujeros
ax3.add_patch(Circle((0, 0), D_HOLE/2, fc='white', ec='black', lw=1.2))
ax3.add_patch(Circle((Lpc, 0), D_HOLE/2, fc='white', ec='black', lw=1.2))
ax3.add_patch(Circle((-BASE_LEN/2, 0), SCREW_DIA/2, fc='white', ec='black', lw=1.2))

# Cotas
draw_dim_line(ax3, (0, 0), (Lpc, 0), 'Lpc = %g' % Lpc, offset_mag=-3.5, fontsize=8)
draw_dim_line(ax3, (-BASE_LEN, -BASE_WID/2 - 0.5), (0, -BASE_WID/2 - 0.5),
               'BASE = %g' % BASE_LEN, offset_mag=-2, fontsize=7, color='#666')
draw_dim_line(ax3, (-BASE_LEN, -BASE_WID/2), (-BASE_LEN, BASE_WID/2),
               'W_b = %g' % BASE_WID, offset_mag=-2, fontsize=7, color='#666')
ax3.annotate('agujero pin %g' % D_HOLE, xy=(Lpc, 0),
             xytext=(Lpc + 1, W_LINK + 1), fontsize=6.5, color='#c62828',
             arrowprops=dict(arrowstyle='->', color='#c62828', lw=0.6))
ax3.annotate('tornillo M%g' % SCREW_DIA, xy=(-BASE_LEN/2, 0),
             xytext=(-BASE_LEN/2, -BASE_WID - 1), fontsize=6.5, color='#c62828',
             arrowprops=dict(arrowstyle='->', color='#c62828', lw=0.6))

ax3.set_xlim(-BASE_LEN - 3, body_len + 4)
ax3.set_ylim(-BASE_WID - 4, W_LINK + 3)
ax3.set_aspect('equal')
ax3.grid(True, ls=':', alpha=0.3)
ax3.set_title('POSTE Lpc (manivela del 4B#3)\n(vista en planta)',
               fontsize=10, fontweight='bold')
ax3.set_xlabel('mm'); ax3.set_ylabel('mm')


# ---- Subplot 2b: Poste Lpc vista lateral ----
ax3b = plt.subplot(3, 2, 4)
ax3b.add_patch(Rectangle((-BASE_LEN, -BASE_THK), BASE_LEN, BASE_THK,
                          fc='#c8e6c9', ec='black', lw=1.5))
ax3b.add_patch(Rectangle((0, 0), body_len, T_LINK,
                          fc='#c8e6c9', ec='black', lw=1.5))
draw_dim_line(ax3b, (0, T_LINK + 0.5), (Lpc, T_LINK + 0.5), 'Lpc = %g' % Lpc,
               offset_mag=2, fontsize=8)
draw_dim_line(ax3b, (-BASE_LEN, -BASE_THK), (-BASE_LEN, 0),
               'T_base = %g' % BASE_THK, offset_mag=-2, fontsize=7, color='#666')
draw_dim_line(ax3b, (body_len, 0), (body_len, T_LINK),
               'T = %g' % T_LINK, offset_mag=2, fontsize=7, color='#666')
ax3b.set_xlim(-BASE_LEN - 4, body_len + 5)
ax3b.set_ylim(-BASE_THK - 3, T_LINK + 5)
ax3b.set_aspect('equal')
ax3b.grid(True, ls=':', alpha=0.3)
ax3b.set_title('Poste Lpc (vista lateral)', fontsize=10)
ax3b.set_xlabel('mm'); ax3b.set_ylabel('mm')


# ---- Subplot 3: Poste Lpd ----
ax4 = plt.subplot(3, 2, 5)
body_len_d = Lpd + D_HOLE / 2 + END_PAD
ax4.add_patch(Rectangle((-BASE_LEN, -BASE_WID/2), BASE_LEN, BASE_WID,
                          fc='#ffcc80', ec='black', lw=1.5))
ax4.add_patch(Rectangle((0, -W_LINK/2), body_len_d, W_LINK,
                          fc='#ffcc80', ec='black', lw=1.5))
ax4.add_patch(Circle((0, 0), D_HOLE/2, fc='white', ec='black', lw=1.2))
ax4.add_patch(Circle((Lpd, 0), D_HOLE/2, fc='white', ec='black', lw=1.2))
ax4.add_patch(Circle((-BASE_LEN/2, 0), SCREW_DIA/2, fc='white', ec='black', lw=1.2))

draw_dim_line(ax4, (0, 0), (Lpd, 0), 'Lpd = %g' % Lpd, offset_mag=-3.5, fontsize=8)
draw_dim_line(ax4, (-BASE_LEN, -BASE_WID/2 - 0.5), (0, -BASE_WID/2 - 0.5),
               'BASE = %g' % BASE_LEN, offset_mag=-2, fontsize=7, color='#666')

ax4.set_xlim(-BASE_LEN - 3, body_len_d + 3)
ax4.set_ylim(-BASE_WID - 3, W_LINK + 3)
ax4.set_aspect('equal')
ax4.grid(True, ls=':', alpha=0.3)
ax4.set_title('POSTE Lpd (balancin del 4B#3)\n(vista en planta)',
               fontsize=10, fontweight='bold')
ax4.set_xlabel('mm'); ax4.set_ylabel('mm')


# ---- Subplot 3b: Cuadro de info / lista de materiales ----
ax5 = plt.subplot(3, 2, 6)
ax5.axis('off')
bom_text = (
    'LISTA DE MATERIALES (BOM):\n'
    '==========================\n\n'
    '  1. Acoplador Lac\n'
    '     - Material: aluminio, acero, PLA o PETG\n'
    '     - Espesor: %g mm\n'
    '     - Cantidad: 1\n\n'
    '  2. Poste Lpc (manivela)\n'
    '     - Rigido sobre falange proximal Fp\n'
    '     - Longitud efectiva: %g mm\n'
    '     - Pivote en IFP, agujero a %g mm del pivote\n'
    '     - Cantidad: 1\n\n'
    '  3. Poste Lpd (balancin)\n'
    '     - Rigido sobre falange distal Fd\n'
    '     - Longitud efectiva: %g mm\n'
    '     - Pivote en IFD, agujero a %g mm del pivote\n'
    '     - Cantidad: 1\n\n'
    '  4. Pasadores (pins)\n'
    '     - Diametro: %g mm\n'
    '     - Cantidad: 2 (CRK3 y ROK3)\n\n'
    '  5. Tornillos para base de montaje\n'
    '     - M%g\n'
    '     - Cantidad: 2 (uno por poste)\n\n'
    'PARAMETROS DE FABRICACION:\n'
    '==========================\n'
    '  Diametro pin (D_PIN)     = %g mm\n'
    '  Diametro agujero (D_HOLE) = %g mm\n'
    '  Holgura radial            = 0.10 mm\n'
    '  Espesor eslabon (T_LINK)  = %g mm\n'
    '  Ancho oreja (W_LINK)      = %g mm\n'
    '  Padding extremo (END_PAD) = %g mm\n'
    '\n'
    'IMPORTANTE:\n'
    '-----------\n'
    'Los planos muestran las piezas en su sistema\n'
    'LOCAL. En el ensamble final, el poste Lpc se\n'
    'rota BETA1_geom = %.1f deg respecto al eje Fp,\n'
    'y el poste Lpd se rota BETA2_geom = %.1f deg\n'
    'respecto al eje Fd. Ver plano general.'
) % (T_LINK, Lpc, Lpc, Lpd, Lpd, D_PIN, SCREW_DIA, D_PIN, D_HOLE, T_LINK,
     W_LINK, END_PAD, beta1_geom, beta2_geom)
ax5.text(0.0, 1.0, bom_text, transform=ax5.transAxes,
         fontsize=8, va='top', ha='left',
         fontfamily='monospace',
         bbox=dict(boxstyle='round,pad=0.4', fc='#fffde7', ec='#999', lw=0.8))

plt.suptitle('Plano de fabricacion - Piezas del 4B#3', fontsize=13,
              fontweight='bold', y=0.995)
plt.tight_layout()
plt.savefig('plano_4B3_piezas.png', dpi=200, bbox_inches='tight',
            facecolor='white')
plt.close()

print("Figura 2 guardada: plano_4B3_piezas.png")
