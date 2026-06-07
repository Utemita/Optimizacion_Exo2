# Analisis de Modificacion de la Articulacion DIP del Exoesqueleto de Mano

> Version 2. Corrige el enfoque de la version 1 (eslabon remoto S2->S3), que
> producia excursiones DIP de ~355 grados y cruce de trayectorias.

## 1. Planteamiento del Problema

La articulacion interfalangica distal (DIP) del exoesqueleto de mano no posee
movimiento independiente. En el modelo cinematico original (`CinematicaExoFinal.m`)
el angulo de la falange distal se calcula como un offset constante respecto a la
falange medial:

```matlab
THETAfd(j) = THETAfm(j) + THETAauxfd;   % THETAauxfd = 38.78 grados (constante)
```

**Consecuencia:** el angulo relativo DIP-PIP es siempre 38.78 grados. La falange
distal queda rigidamente unida a la medial y, en el modelo CAD, no presenta
articulacion propia. Las tres falanges deben flexionarse unas respecto a otras.

## 2. Analisis de Causa Raiz

La cadena cinematica transmite el movimiento desde un unico motor hasta las
falanges:

| Etapa | Mecanismo | Salida | DOF independiente |
|-------|-----------|--------|-------------------|
| 1 | Par de engranes | theta_1, theta_2 | relacion fija reng |
| 2 | 5 barras #1 | Punto P | si |
| 3 | 4 barras #1 | theta_fp (falange proximal) | si |
| 4 | 5 barras #2 | Punto P2 | si |
| 5 | 4 barras #2 | theta_fm (falange medial) | si |
| 6 | **Offset constante** | **theta_fd = theta_fm + 38.78** | **NO** |

Los mecanismos 1-5 dan movimiento variable a las falanges proximal y medial,
pero la falange distal solo hereda el angulo de la medial mas una constante.
No existe ningun mecanismo fisico que transmita movimiento diferencial al DIP.

## 3. Solucion Propuesta (v2) - Mecanismo de 4 Barras que Cruza la Articulacion IFP

### 3.1 Concepto

Se agrega un mecanismo de 4 barras de **acoplamiento** que **cruza la articulacion
IFP (PIP)**, utilizando la **falange medial como bancada (ground)**. Su longitud de
bancada es `fm`, que es conocida y fija.

- **Bancada (d3 = fm):** la propia falange medial, de la articulacion IFP (O1) a la
  articulacion IFD (O2).
- **Manivela de entrada (Lpc):** poste rigido montado sobre la falange **proximal**,
  con pivote en la articulacion IFP (O1). Su orientacion absoluta es
  `theta_fp + BETA1`, es decir, la impone directamente la falange proximal.
- **Balancin de salida (Lpd):** poste rigido montado sobre la falange **distal**, con
  pivote en la articulacion IFD (O2). Su orientacion determina `theta_fd`.
- **Acoplador (Lac):** barra que une el extremo de la manivela con el del balancin.

### 3.2 Principio de Funcionamiento

Conforme la articulacion IFP (PIP) se flexiona, el angulo relativo entre la falange
proximal y la medial cambia (en los datos del mecanismo varia unos 73 grados). Como
la manivela esta solidaria a la proximal y la bancada es la medial, ese movimiento
relativo acciona el mecanismo de 4 barras, que a traves del acoplador hace girar el
balancin (solidario a la distal) y produce una flexion DIP **variable**, monotona y
acotada. Todo con un **unico motor** (no se anaden actuadores).

### 3.3 Diagrama de Topologia

```
        Falange PROXIMAL                 Falange MEDIAL = BANCADA              Falange DISTAL
        (impone theta_fp)                (longitud fm, de O1 a O2)            (define theta_fd)

                 *  punta manivela                       * punta balancin
                / \                                     / \
               /   \  Lpc                              /   \  Lpd
              /     \ (poste prox.)        Lac        /     \  (poste dist.)
             /       \ - - - - - - - - - - - - - - - /       \
            /         *  - - - acoplador - - -  *             \
           /        O1 = IFP                    O2 = IFD        \
   =======*=================================================*===========
          |<--------------------- fm (bancada) ------------->|
        BETA1 = montaje manivela vs proximal      BETA2 = montaje balancin vs distal
```

### 3.4 Ventajas frente a la version 1

- No usa el "eslabon remoto S2->S3" con longitud de manivela fija desligada de la
  distancia real (esa inconsistencia provocaba las vueltas de 355 grados).
- La bancada `fm` es una longitud real y fija: no introduce parametros geometricos
  redundantes.
- La salida se resuelve con una sola rama de ensamble y desenrollado de continuidad,
  evitando saltos por cambio de rama o por el corte del arcotangente.

## 4. Ecuaciones Cinematicas del Tercer Mecanismo

### 4.1 Marco de referencia local de la falange medial

- Origen en O1 = IFP. Eje X local a lo largo de la falange medial (orientacion
  `theta_fm`), eje Y perpendicular.
- O1_local = (0, 0); O2_local = (fm, 0).

### 4.2 Manivela de entrada

El angulo de la manivela en el marco local de la medial es:

```
alpha1 = (theta_fp + BETA1) - theta_fm
```

Punta de la manivela (marco local):

```
A = ( Lpc*cos(alpha1) , Lpc*sin(alpha1) )
```

### 4.3 Resolucion del acoplador (1 ecuacion, 1 incognita)

El balancin tiene punta `B = O2 + Lpd*(cos(alpha2), sin(alpha2))`, con O2 = (fm, 0).
La restriccion del acoplador es `|A - B| = Lac`. Definiendo el vector de O2 a A:

```
Px = A_x - fm
Py = A_y
R  = sqrt(Px^2 + Py^2)
```

al desarrollar `|A - B|^2 = Lac^2` se obtiene la ecuacion lineal en seno-coseno:

```
Px*cos(alpha2) + Py*sin(alpha2) = K
con  K = (Px^2 + Py^2 + Lpd^2 - Lac^2) / (2*Lpd)
```

cuya solucion es:

```
phi    = atan2(Py, Px)
alpha2 = phi - acos(K / R)        % rama de ensamble (configuracion abierta)
```

**Condicion de ensamble:** `|K| <= R`. Si no se cumple, el mecanismo no se puede
ensamblar en esa posicion y se usa el offset constante como respaldo (fallback).

Para garantizar continuidad se aplica un desenrollado del angulo respecto al paso
anterior (se suma o resta 2*pi si el salto supera pi).

### 4.4 Angulo de salida de la falange distal

El balancin esta solidario a la falange distal con un offset de montaje BETA2, de
modo que `alpha2 = (theta_fd + BETA2) - theta_fm`. Despejando:

```
theta_fd = theta_fm + alpha2 - BETA2
```

`theta_fd` ya NO es constante respecto a `theta_fm`: varia conforme se flexiona el PIP.

## 5. Tabla de Parametros del Tercer Mecanismo

| Parametro | Simbolo | Significado | Valor inicial | Unidad |
|-----------|---------|-------------|---------------|--------|
| Manivela  | Lpc   | Poste sobre la falange proximal (pivote IFP) | 8     | mm |
| Balancin  | Lpd   | Poste sobre la falange distal  (pivote IFD)  | 18    | mm |
| Acoplador | Lac   | Barra que une manivela y balancin            | 8.86  | mm |
| Montaje 1 | BETA1 | Angulo manivela respecto a falange proximal  | 40    | grados |
| Montaje 2 | BETA2 | Angulo balancin respecto a falange distal    | 110   | grados |
| Bancada   | fm    | Falange medial (NO es parametro libre)       | 26    | mm |

La bancada es `fm` (longitud real de la falange medial), por lo que no es un
parametro de diseno independiente.

### 5.1 Montaje fisico DORSAL del 4B#3 y equivalencia palmar-dorsal

Por restriccion fisica el mecanismo del DIP se construye del lado **dorsal** del
dedo (lado opuesto a la palma), porque las falanges apoyan el lado palmar contra
los objetos manipulados y el exoesqueleto no puede invadir esa zona. Las tres
articulaciones (MCF, IFP/PIP, IFD/DIP) quedan asi sostenidas por mecanismos
dorsales (5B#1 sobre MCF, 5B#2 + 4B#1 sobre la proximal-IFP, 4B#2 sobre la medial,
y el nuevo 4B#3 sobre la distal-IFD).

**Propiedad clave (simetria del 4 barras):** el 4 barras plano es invariante bajo
reflexion respecto a su bancada. Para el 4B#3, cuya bancada es la falange medial
(eje IFP-IFD), un mecanismo "palmar" y su reflejo "dorsal" tienen las mismas
longitudes (Lpc, Lpd, Lac), los mismos offsets de montaje (BETA1, BETA2) en valor
absoluto, y producen **identico angulo de salida** `theta_fd`. La unica diferencia
es donde quedan fisicamente las barras intermedias (manivela, acoplador, balancin):
arriba (dorsal) o abajo (palmar) de la falange medial.

**Consecuencia para el codigo:** las ecuaciones cinematicas de la Seccion 4 son
validas tal cual, sin cambio alguno. La cinematica numerica
(`CinematicaExoModificada.m`, `exo_18.py`) calcula `theta_fd` correctamente con la
convencion estandar. Solo el **diagrama de eslabones** refleja las posiciones de la
manivela y el balancin respecto a la recta IFP-IFD (la bancada) para mostrar el
mecanismo donde se construira realmente en el CAD.

**Implementacion del reflejo en `diagrama_mecanismo.py`:** dado un punto P palmar,
su reflejo dorsal respecto a la recta que pasa por A en direccion unitaria u_hat es

```
v        = P - A
parallel = (v . u_hat) * u_hat
perp     = v - parallel
P_dorsal = A + parallel - perp        # invierte la componente perpendicular
```

Aplicado a los puntos `CRK3` (extremo de la manivela) y `ROK3` (extremo del
balancin) con `A = IFP` y `u_hat = (IFD - IFP) / fm`, los lleva al lado dorsal sin
modificar la cinematica de `theta_fd`.

## 6. Criterios de Verificacion (resultados con los valores iniciales)

Simulando la cadena completa con los parametros nominales del dedo indice:

| Criterio | Objetivo | Resultado obtenido |
|----------|----------|--------------------|
| Excursion DIP | 20-40 grados | **29.9 grados** (de 32.1 a 62.0) |
| Monotonia | sin inversiones | **monotono creciente** |
| Sin cruce de trayectorias IFP/IFD/punta | no se cruzan | **no se cruzan** |
| Ensamble en todo el rango | margen `R - |K| > 0` | **margen minimo 0.086 (>0)** |
| Acoplamiento DIP/PIP | aprox 0.3-0.7 | ~0.4 (29.9/73.4) |

El angulo DIP relativo deja de ser la linea constante de 38.78 grados y pasa a ser
una curva creciente, lo que reproduce la flexion progresiva de la falange distal
durante el cierre del dedo.

## 7. Codigo MATLAB (fragmentos clave)

### Definicion de parametros (al inicio del script):

```matlab
Lpc   = 8;     % Manivela: poste sobre la falange PROXIMAL (pivote en IFP) [mm]
Lpd   = 18;    % Balancin: poste sobre la falange DISTAL  (pivote en IFD) [mm]
Lac   = 8.86;  % Acoplador que une las puntas de manivela y balancin [mm]
BETA1 = 40;    % Angulo de montaje de la manivela respecto a la prox. [grados]
BETA2 = 110;   % Angulo de montaje del balancin respecto a la distal [grados]
```

### Dentro del bucle (despues de calcular pxIFD, pyIFD):

```matlab
% Manivela (en la proximal) expresada en el marco local de la medial
alpha1_3 = thetafp(j) + deg2rad(BETA1) - thetafm(j);
Ax3 = Lpc*cos(alpha1_3);
Ay3 = Lpc*sin(alpha1_3);

% Restriccion del acoplador: |A - B| = Lac, con O2 = (fm, 0)
Px3 = Ax3 - fm;
Py3 = Ay3;
R3  = sqrt(Px3^2 + Py3^2);
K3  = (Px3^2 + Py3^2 + Lpd^2 - Lac^2) / (2*Lpd);

if abs(K3) > R3
    warning('Tercer mecanismo no ensambla, paso %d', j);
    THETAfd(j) = THETAfm(j) + THETAauxfd;          % fallback
else
    phi3 = atan2(Py3, Px3);
    alpha2_3 = phi3 - acos(K3 / R3);               % rama de ensamble
    if exist('alpha2_3_prev', 'var')               % continuidad
        while (alpha2_3 - alpha2_3_prev) >  pi, alpha2_3 = alpha2_3 - 2*pi; end
        while (alpha2_3 - alpha2_3_prev) < -pi, alpha2_3 = alpha2_3 + 2*pi; end
    end
    alpha2_3_prev = alpha2_3;
    THETAfd(j) = THETAfm(j) + rad2deg(alpha2_3 - deg2rad(BETA2));
end

thetafd(j) = deg2rad(THETAfd(j));
pxPF(j) = fd*cos(thetafd(j)) + pxIFD(j);
pyPF(j) = fd*sin(thetafd(j)) + pyIFD(j);
```

La implementacion completa esta en `CinematicaExoModificada.m`.

## 8. Impacto en la Optimizacion Python (exo_18.py)

### 8.1 Nuevos parametros de optimizacion

Se agregan 5 parametros al vector de diseno (los 17 actuales pasan a 22):

| Indice | Parametro | Valor inicial | Limite inf. | Limite sup. |
|--------|-----------|---------------|-------------|-------------|
| 17 | Lpc   | 0.008 m | 0.005 | 0.025 |
| 18 | Lpd   | 0.018 m | 0.008 | 0.030 |
| 19 | Lac   | 0.0089 m | 0.004 | 0.030 |
| 20 | BETA1 | 0.698 rad (40 grados) | 0.0  | 2.97 (170 grados) |
| 21 | BETA2 | 1.920 rad (110 grados) | 0.0 | 2.97 (170 grados) |

### 8.2 Reemplazo en `run_kinematics`

Sustituir `theta_fd = theta_fm + theta_aux_fd` por:

```python
# Tercer mecanismo de 4 barras (bancada = falange medial, longitud FM_REAL)
alpha1_3 = theta_fp + beta1 - theta_fm        # manivela en marco local medial
Ax3 = Lpc * np.cos(alpha1_3)
Ay3 = Lpc * np.sin(alpha1_3)
Px3 = Ax3 - FM_REAL
Py3 = Ay3
R3  = np.hypot(Px3, Py3)
K3  = (Px3**2 + Py3**2 + Lpd**2 - Lac**2) / (2.0 * Lpd)
if abs(K3) > R3:
    return None                                # no ensambla -> penalizar
phi3 = np.arctan2(Py3, Px3)
alpha2_3 = phi3 - np.arccos(K3 / R3)           # rama de ensamble
# (desenrollar respecto al paso previo para continuidad)
theta_fd = theta_fm + (alpha2_3 - beta2)
```

### 8.3 Consideraciones

1. **Penalizacion por no ensamble:** si `|K3| > R3` en algun paso, descartar el
   conjunto de parametros (devolver None / fitness alto).
2. **Continuidad:** desenrollar `alpha2_3` respecto al valor del paso anterior para
   evitar saltos por el corte del arcocoseno/arcotangente.
3. **Monotonia:** la penalizacion anti-gancho existente puede extenderse al angulo
   DIP relativo (`theta_fd - theta_fm`) para forzar flexion progresiva.
4. **Inicializacion:** usar los valores de la Seccion 5, que ya ensamblan en todo el
   rango y dan ~30 grados de excursion DIP monotona.
