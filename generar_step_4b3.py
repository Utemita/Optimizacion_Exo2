"""
Genera archivos STEP de las piezas nuevas del tercer mecanismo de 4 barras (4B#3).

Piezas:
  1. acoplador_Lac.step      - barra plana con dos agujeros (8.86 mm entre centros)
  2. poste_Lpc.step          - oreja para montar sobre la falange proximal
  3. poste_Lpd.step          - oreja para montar sobre la falange distal
  4. pasador.step            - pin de articulacion (referencia)

Convenciones:
  - Plano del mecanismo: XY (Z = espesor).
  - Origen de cada pieza coincide con el punto pivote principal.
  - Para los postes, el "frente" (donde sale el agujero del pin) esta en +X.
  - Las dimensiones siguen los parametros del 4B#3.

Tolerancias y convenciones de manufactura:
  - Diametro nominal del pasador: D_PIN = 1.5 mm (default).
  - Holgura agujero-pin: HOLE_CLEARANCE = 0.10 mm radial.
  - Espesor de las orejas y del acoplador: T_LINK = 1.5 mm.
  - Anchura de las orejas: W_LINK = 4 mm.
"""
import cadquery as cq
import os
import math

# ===================== PARAMETROS GLOBALES (mm) =====================
# Cinematicos (del modelo)
LPC = 8.0          # longitud manivela (poste sobre Fp)
LPD = 18.0         # longitud balancin (poste sobre Fd)
LAC = 8.86         # longitud acoplador

# Constructivos (configurables)
D_PIN = 1.5                # diametro nominal pin
HOLE_CLEARANCE = 0.10      # holgura radial agujero-pin
D_HOLE = D_PIN + 2 * HOLE_CLEARANCE  # diametro de los agujeros
T_LINK = 1.5               # espesor de eslabones/orejas
W_LINK = 4.0               # ancho nominal de eslabones/orejas
END_PAD = 1.5              # material desde el borde del agujero al final del eslabon
ROOT_PAD = 1.5             # material en la base del poste

# Base de montaje (donde el poste se atornilla a la falange)
BASE_LEN = 6.0             # longitud de la base de montaje a lo largo de la falange
BASE_WID = 4.0             # ancho de la base de montaje (perpendicular al eje del dedo)
BASE_THK = 1.0             # grueso de la base
SCREW_DIA = 1.5            # diametro de los tornillos de montaje (M1.5 o similar)

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
STEP_DIR = os.path.join(OUTPUT_DIR, 'step_4B3')
os.makedirs(STEP_DIR, exist_ok=True)


# ===================== ACOPLADOR Lac =====================
def make_acoplador():
    """Barra plana con dos agujeros pasantes a una distancia LAC entre centros.
    Eje de la barra a lo largo de X, espesor en Z, agujeros en Z.
    """
    L_total = LAC + 2 * (D_HOLE / 2 + END_PAD)
    body = (cq.Workplane('XY')
            .box(L_total, W_LINK, T_LINK, centered=(True, True, True)))
    body = (body.faces('>Z')
            .workplane()
            .pushPoints([(-LAC / 2, 0), (LAC / 2, 0)])
            .hole(D_HOLE))

    return body


# ===================== POSTE Lpc (manivela) =====================
def make_poste(L_poste, mount_angle_deg, name='poste'):
    """Poste tipo oreja con base de montaje.

    El poste arranca del pivote (origen) y sale a longitud L_poste con un agujero
    en el extremo donde se conecta el acoplador.

    Para simplicidad, generamos la oreja en su sistema local: el poste sale a lo
    largo de +X. El usuario rotara la pieza al ensamblarla en CAD para coincidir
    con BETA1_geom o BETA2_geom segun corresponda.

    Args:
        L_poste: longitud del poste (Lpc o Lpd) [mm]
        mount_angle_deg: angulo de montaje (referencia, no se aplica aqui;
                         el ensamblaje final se hace en el CAD del usuario)
        name: nombre identificador (no se usa en la geometria)
    """
    # Cuerpo de la oreja: rectangulo desde el pivote hasta el extremo del poste
    body_len = L_poste + D_HOLE / 2 + END_PAD
    body = (cq.Workplane('XY')
            .box(body_len, W_LINK, T_LINK, centered=(False, True, True)))
    # Agujero en el pivote (en X=0)
    body = (body.faces('>Z')
            .workplane()
            .pushPoints([(0, 0)])
            .hole(D_HOLE))
    # Agujero en el extremo (X = L_poste)
    body = (body.faces('>Z')
            .workplane(centerOption='CenterOfMass')
            .pushPoints([(L_poste - body_len / 2, 0)])
            .hole(D_HOLE))

    # Base de montaje: una pestana plana en la cara opuesta del pivote (X negativo)
    # que se atornilla a la falange. La base sale en -X desde el pivote.
    base = (cq.Workplane('XY')
            .center(-BASE_LEN / 2, 0)
            .box(BASE_LEN, BASE_WID, BASE_THK, centered=(True, True, True))
            .translate((0, 0, -T_LINK / 2 - BASE_THK / 2)))
    # Agujero pasante para tornillo en la base (un solo tornillo central)
    base = (base.faces('>Z')
            .workplane()
            .pushPoints([(0, 0)])
            .hole(SCREW_DIA))

    # Combinar
    assembly = body.union(base)
    return assembly


def make_poste_lpc():
    return make_poste(LPC, 0, 'Lpc')


def make_poste_lpd():
    return make_poste(LPD, 0, 'Lpd')


# ===================== PASADOR =====================
def make_pasador():
    """Pin cilindrico de diametro D_PIN, longitud suficiente para atravesar
    dos orejas de espesor T_LINK con un poco de holgura."""
    L_pin = 2 * T_LINK + 1.0  # 1 mm de holgura total
    pin = (cq.Workplane('XY')
           .circle(D_PIN / 2)
           .extrude(L_pin)
           .faces('>Z')
           .workplane()
           .circle(D_PIN / 2 + 0.4)  # cabeza tipo "T" (opcional, pequena)
           .extrude(0.4))
    return pin


# ===================== EXPORTAR =====================
print("Generando archivos STEP...")
print(f"  Directorio de salida: {STEP_DIR}")
print()

print("[1/4] Acoplador Lac (%.2f mm entre centros, %.2f de espesor)..." % (LAC, T_LINK))
acoplador = make_acoplador()
acoplador.val().exportStep(os.path.join(STEP_DIR, 'acoplador_Lac.step'))
print(f"      Guardado: acoplador_Lac.step")

print(f"\n[2/4] Poste manivela Lpc (longitud {LPC} mm)...")
lpc_poste = make_poste_lpc()
lpc_poste.val().exportStep(os.path.join(STEP_DIR, 'poste_Lpc.step'))
print(f"      Guardado: poste_Lpc.step")

print(f"\n[3/4] Poste balancin Lpd (longitud {LPD} mm)...")
lpd_poste = make_poste_lpd()
lpd_poste.val().exportStep(os.path.join(STEP_DIR, 'poste_Lpd.step'))
print(f"      Guardado: poste_Lpd.step")

print(f"\n[4/4] Pasador (diametro {D_PIN} mm)...")
pin = make_pasador()
pin.val().exportStep(os.path.join(STEP_DIR, 'pasador.step'))
print(f"      Guardado: pasador.step")

# Tambien exportamos un STL de cada pieza (formato mas compatible)
print("\nExportando STL adicionales (mas compatibles)...")
acoplador.val().exportStl(os.path.join(STEP_DIR, 'acoplador_Lac.stl'))
lpc_poste.val().exportStl(os.path.join(STEP_DIR, 'poste_Lpc.stl'))
lpd_poste.val().exportStl(os.path.join(STEP_DIR, 'poste_Lpd.stl'))
pin.val().exportStl(os.path.join(STEP_DIR, 'pasador.stl'))
print("  Listo.")

print("\n=== Resumen de piezas ===")
print(f"  Acoplador Lac:    {LAC:.2f} mm entre agujeros, {W_LINK} mm ancho, {T_LINK} mm espesor")
print(f"  Poste Lpc:        {LPC} mm efectivo, base de montaje {BASE_LEN}x{BASE_WID}x{BASE_THK} mm")
print(f"  Poste Lpd:        {LPD} mm efectivo, base de montaje {BASE_LEN}x{BASE_WID}x{BASE_THK} mm")
print(f"  Pasador (pin):    diametro {D_PIN} mm, longitud {2*T_LINK + 1.0} mm")
print(f"  Diametro agujeros: {D_HOLE} mm (holgura radial {HOLE_CLEARANCE} mm)")
