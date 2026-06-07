# Guia de Ensamble CAD - Exoesqueleto de Dedo Indice

Este documento describe paso a paso como construir el modelo CAD del exoesqueleto
a partir del diagrama de eslabones generado por `diagrama_mecanismo.py`. Seguir
esta guia garantiza que el estudio de movimiento en CAD reproduzca las mismas
trayectorias que `CinematicaExoModificada.m`.

---

## 1. Sistema de referencia

| Elemento | Posicion (mm) |
|----------|---------------|
| Origen global | Punto A = (-9, 0) |
| Punto B | (9, 0) |
| Art. MCF | (-9, -20) |

El eje X apunta hacia la derecha (dorsal-distal), el eje Y hacia arriba.
La bancada es un triangulo rigido formado por A, B y MCF.

---

## 2. Tabla completa de parametros

| Parametro | Simbolo | Valor | Unidad | Descripcion |
|-----------|---------|-------|--------|-------------|
| Bancada1 | 2*r3 | 18 | mm | Distancia entre pivotes A y B |
| Bancada2 | d | 20 | mm | Distancia de A a MCF (vertical) |
| Link1 | L1 (= r4) | 35 | mm | Manivela del 5B#1 desde B |
| Link2 | L2 (= r5) | 49 | mm | Eslabon J2-P del 5B#1 |
| Link3 | L3 (= r2) | 25 | mm | Eslabon M4-P del 5B#1 |
| Link4 | L4 (= r1 = a) | 20 | mm | Manivela del 5B#1 desde A / manivela 4B#1 |
| Link5 | L5 (= b) | 25 | mm | Eslabon M4-S1 del 4B#1 |
| Link6 | L6 (= r5m2) | 55 | mm | Eslabon P-P2 del 5B#2 |
| Link7 | L7 (= r2m2 = a2) | 35 | mm | Eslabon S2-P2 del 5B#2 / manivela 4B#2 |
| Link8 | L8 (= b2) | 52 | mm | Acoplador P2-P3 del 4B#2 |
| c2 | c2 | 46.01 | mm | Balancin IFP-P3 del 4B#2 |
| hsp | hsp | 17 | mm | Altura de los soportes (perpendicular a falange) |
| dsp | dsp | 18 | mm | Distancia del soporte a la articulacion MCF/IFP |
| fp | Fp | 49 | mm | Longitud de la falange proximal |
| fm | Fm | 26 | mm | Longitud de la falange medial |
| fd | Fd | 24 | mm | Longitud de la falange distal |
| THETA1_ini | - | 109 | deg | Angulo inicial de la manivela 4B#1 |
| THETA14B | - | 90 | deg | Angulo de la bancada del 4B#1 |
| THETAauxfm | - | 51.39 | deg | Offset angular entre c2 y Fm |
| reng | - | 2 | - | Relacion de engranaje |
| Lpc | - | 8 | mm | Manivela del 4B#3 (montada en Fp, pivote en IFP) |
| Lpd | - | 18 | mm | Balancin del 4B#3 (montado en Fd, pivote en IFD) |
| Lac | - | 8.86 | mm | Acoplador del 4B#3 |
| BETA1_geom | - | -11.4 | deg | **Angulo geometrico** de Lpc respecto al eje Fp (medido en el dibujo dorsal) |
| BETA2_geom | - | 185.8 | deg | **Angulo geometrico** de Lpd respecto al eje Fd (medido en el dibujo dorsal) |

> **Nota sobre los angulos BETA1 y BETA2:** En el codigo `CinematicaExoModificada.m`
> aparecen como `BETA1=40` y `BETA2=110`. Estos son los valores que entran en la
> formulacion analitica (lado palmar). En el diagrama dorsal real (donde tu vas
> a construir el mecanismo) los angulos que mide el transportador son
> `BETA1_geom = -11.4 deg` y `BETA2_geom = 185.8 deg`. Las longitudes Lpc, Lpd y
> Lac son las mismas en ambas representaciones, y las trayectorias son identicas.

---

## 3. Topologia del mecanismo (cadena cinematica)

```
MOTOR (en B)
  |
  +-- Engranajes (reng=2) --> Eje en A
  |
  +-- [5-BARRAS #1]  (entradas: th2 en B, th1 en A)
  |     Eslabones: L1(B-J2), L2(J2-P), L3(M4-P), L4(A-M4)
  |     Salida: posicion del punto P
  |
  +-- [4-BARRAS #1]  (entrada: th1 = angulo de L4)
  |     Eslabones: a=L4(manivela), b=L5(M4-S1), c=sqrt(hsp^2+dsp^2)(S1-MCF), d=B2(bancada)
  |     Salida: angulo theta4a --> orientacion falange proximal (THETAfp)
  |
  +-- [5-BARRAS #2]  (entradas: angulos de S1-S2 y M4-P, en ref. rotada)
  |     Eslabones: fp-2dsp(S1-S2), L7(S2-P2), L3/2(mitad5B), L3(otra mitad), L6(P-P2)
  |     Salida: posicion del punto P2
  |
  +-- [4-BARRAS #2]  (entrada: angulo de L7 desde S2)
  |     Eslabones: a2=L7(manivela S2-P2), b2=L8(P2-P3), c2(IFP-P3), d2=sqrt(hsp^2+dsp^2)(IFP-S2)
  |     Salida: angulo de c2 --> THETAfm = THETA4am2 + THETAauxfm
  |     >> Conexion rigida P3-Fm via L9 + poste hsp sobre Fm <<
  |
  +-- [4-BARRAS #3]  (entrada: rotacion relativa Fp-Fm, o sea flexion PIP)
        Bancada: falange medial Fm (IFP a IFD, longitud 26 mm)
        Manivela: Lpc=8 (rigida a Fp, pivote en IFP, BETA1_geom = -11.4 deg)
        Acoplador: Lac=8.86 (eslabon flotante CRK3-ROK3)
        Balancin: Lpd=18 (rigida a Fd, pivote en IFD, BETA2_geom = 185.8 deg)
        Salida: angulo de la falange distal THETAfd
```

---

## 4. Instrucciones de ensamble paso a paso

### 4.1 Marco fijo (Ground)

1. Crear un cuerpo rigido con tres puntos de articulacion: A, B, MCF.
2. A esta en (-9, 0), B en (9, 0), MCF en (-9, -20).
3. Fijar este cuerpo al suelo (restriccion de posicion completa).
4. En B colocar un motor rotativo (entrada THETA2).
5. En A colocar un engranaje que gira a THETA2/2 + 109 deg (relacion 2:1).

### 4.2 Primer mecanismo de 5 barras

1. Crear eslabon L1 (35 mm): pivote revoluta en B, extremo libre = J2.
2. Crear eslabon L4 (20 mm): pivote revoluta en A, extremo libre = M4.
3. L1 gira con THETA2 (motor). L4 gira con THETA1 = THETA2/2 + 109 deg.
4. Crear eslabon L2 (49 mm): conectar J2 con P (revoluta en ambos extremos).
5. Crear eslabon L3 (25 mm): conectar M4 con P (revoluta en ambos extremos).
6. El punto P queda determinado como la interseccion de los circulos desde J2 y M4. Usar la solucion con signo (+) del discriminante (configuracion positiva).

### 4.3 Primer mecanismo de 4 barras (soporte de la falange proximal)

1. Crear eslabon L5 (25 mm): pivote en M4, extremo libre = S1.
2. Crear eslabon c = sqrt(17^2 + 18^2) = 24.76 mm: pivote en MCF, extremo libre = S1.
3. Juntar el extremo de L5 con el extremo de c en S1 (revoluta).
4. La manivela del 4B#1 es L4 (el mismo eslabon que en el 5B#1).
5. La bancada del 4B#1 es la vertical A-MCF (= Bancada2 = 20 mm, angulo 90 deg desde horizontal).
6. Resultado: S1 es el pivote inferior del soporte sobre la falange proximal.

### 4.4 Falange proximal

1. La falange proximal es un cuerpo rigido de longitud fp = 49 mm.
2. Su pivote esta en MCF (articulacion MCP).
3. Su angulo es: THETAfp = THETA4a + atan2(hsp, dsp) = THETA4a + 43.34 deg.
4. Los puntos S1 y S2 estan sobre la falange proximal:
   - S1 esta a distancia c = 24.76 mm de MCF con angulo theta4a.
   - S2 esta a distancia rs2 = sqrt(17^2 + 31^2) = 35.36 mm de MCF, con angulo THETAfp - atan2(17, 31).
5. La distancia entre S1 y S2 a lo largo de la falange es fp - 2*dsp = 13 mm.
6. El extremo distal de la falange proximal es IFP (articulacion PIP).

### 4.5 Segundo mecanismo de 5 barras

Este mecanismo opera en un sistema de referencia rotado cuyo centro esta en el punto medio entre S1 y M4.

1. Crear eslabon L7 (35 mm): pivote en S2, extremo libre = P2.
2. Crear eslabon L6 (55 mm): pivote en P, extremo libre = P2.
3. Conectar ambos extremos en P2 (revoluta).
4. Las "entradas" del 5B#2 son los angulos de S1-S2 y M4-P, transformados al sistema rotado.

### 4.6 Segundo mecanismo de 4 barras (driver de la falange medial)

1. Bancada del 4B#2: segmento IFP-S2 (longitud d2 = sqrt(17^2+18^2) = 24.76 mm).
2. Manivela: L7 (= a2 = 35 mm), pivote en S2.
3. Acoplador: L8 (= b2 = 52 mm), conecta P2 con P3.
4. Balancin: c2 (= 46.01 mm), pivote en IFP, extremo libre = P3.
5. Conectar P2-P3 (revoluta) y P3-IFP (revoluta).
6. El angulo de c2 determina la orientacion de la falange medial:
   THETAfm = THETA4am2 + THETAauxfm (= angulo_balancin + 51.39 deg).

### 4.7 Falange medial y su poste hsp

1. Cuerpo rigido de longitud fm = 26 mm, pivote en IFP.
2. Su angulo es THETAfm. Su extremo distal es IFD.
3. **Poste hsp sobre Fm**: extrusion rigida de longitud hsp = 17 mm, perpendicular
   a Fm en su lado dorsal, anclada a la mitad de Fm (a 13 mm de IFP).
4. **Eslabon L9**: barra rigida que conecta el extremo del poste hsp con P3.
   Esta barra fuerza el offset angular `THETAauxfm = 51.39 deg` entre c2 y Fm.
   En CAD, L9 puede modelarse como una restriccion rigida en lugar de una pieza
   real (mismo efecto: P3 y el extremo del poste hsp tienen distancia constante).

### 4.8 Tercer mecanismo de 4 barras (driver DIP) - LADO DORSAL

Este es el mecanismo nuevo, que reemplaza al offset rigido entre Fd y Fm. Se
construye del lado **DORSAL** del dedo (lado opuesto a la palma).

1. **Bancada**: la propia falange medial (fm = 26 mm, de IFP a IFD).
2. **Manivela Lpc** (8 mm): rigidamente unida a la falange PROXIMAL (no a Fm).
   - Pivote en la articulacion IFP.
   - Anclaje a Fp con `BETA1_geom = -11.4 deg` (casi alineada con Fp, ligeramente
     hacia dorsal). En el dibujo, Lpc sale de IFP "hacia atras" (hacia MCF) con
     una pequena inclinacion dorsal.
   - Su extremo es CRK3.
3. **Acoplador Lac** (8.86 mm): eslabon plano con dos agujeros pasantes.
   - Pivote (revoluta) en CRK3 con la manivela.
   - Pivote (revoluta) en ROK3 con el balancin.
4. **Balancin Lpd** (18 mm): rigidamente unida a la falange DISTAL (no a Fm).
   - Pivote en la articulacion IFD.
   - Anclaje a Fd con `BETA2_geom = 185.8 deg` (casi opuesta al eje Fd, dorsal).
     Es decir, Lpd sale de IFD apuntando "hacia IFP" (proximalmente) con una
     pequena inclinacion dorsal.
   - Su extremo es ROK3.

**Como funciona:** Cuando la articulacion PIP se flexiona (la falange medial rota
respecto a la proximal), la manivela Lpc -- que esta rigida a la proximal -- cambia
su angulo relativo respecto a la bancada (la medial). Este cambio angular acciona
el mecanismo de 4 barras y produce una rotacion del balancin Lpd, que a su vez
mueve la falange distal.

**Resultado:** la articulacion DIP se mueve con una excursion de ~30 grados
(32 a 62 deg relativo a la medial) de forma monotona durante el barrido completo
de la manivela de entrada.

### 4.9 Falange distal

1. Cuerpo rigido de longitud fd = 24 mm.
2. Pivote en IFD (articulacion DIP).
3. Su angulo es THETAfd, determinado por el 4B#3.
4. Su punta es el extremo del dedo (TIP).
5. **Poste hsp sobre Fd**: en el mecanismo original llevaba L10. En el diseno
   nuevo este poste queda VESTIGIAL (no se conecta a nada). Puede omitirse en
   el CAD o dejarse como referencia geometrica.

---

## 5. Restricciones para el estudio de movimiento

| Tipo | Ubicacion | Descripcion |
|------|-----------|-------------|
| Motor rotativo | Pivote B | Entrada THETA2 (0 a 132 deg para cierre completo) |
| Engranaje | A respecto a B | THETA1 = THETA2/reng + 109 deg |
| Revoluta | Cada union de eslabones | Permite rotacion relativa |
| Rigida | Lpc con Fp | La manivela del 4B#3 es parte del cuerpo de la falange proximal |
| Rigida | Lpd con Fd | El balancin del 4B#3 es parte del cuerpo de la falange distal |
| Rigida | Poste hsp con Fm | El poste sobre Fm es parte del cuerpo de la falange medial |
| Rigida | Poste hsp con Fd | El poste sobre Fd es parte del cuerpo de la falange distal |
| Rigida | L9 con poste hsp(Fm) y P3 | Materializa el offset THETAauxfm |
| Fijo | Marco A-B-MCF | Completamente fijo al suelo |

---

## 6. Verificacion del ensamble

Despues de construir el modelo CAD y correr el estudio de movimiento:

1. **Trayectoria IFP**: debe coincidir con la curva proximal del MATLAB.
2. **Trayectoria IFD**: debe coincidir con la curva medial del MATLAB.
3. **Trayectoria TIP**: debe coincidir con la curva distal del MATLAB.
4. **Angulo DIP relativo (THETAfd - THETAfm)**: debe ir de ~32 a ~62 grados.
5. **Monotonia**: el angulo DIP relativo debe ser estrictamente creciente (sin inversiones).

Para generar las trayectorias de referencia en formato CSV, ejecutar:
```bash
python3 generar_trayectorias_referencia.py
```
Esto produce `trayectorias_referencia.csv` con las columnas:
THETA2, pxIFP, pyIFP, pxIFD, pyIFD, pxTIP, pyTIP, DIP_relativo

---

## 7. Notas sobre el montaje fisico

- **Diagrama generado** (`diagrama_mecanismo_completo.png`): muestra las
  posiciones EXACTAS calculadas por las ecuaciones de MATLAB, con CRK3 y ROK3
  reflejados al lado dorsal.
- **Plano general cotado** (`plano_4B3_general.png`): vista lateral del 4B#3
  con `BETA1_geom`, `BETA2_geom`, Lpc, Lac, Lpd cotados.
- **Plano de piezas** (`plano_4B3_piezas.png`): planos de fabricacion
  individuales del acoplador y los postes con todas las dimensiones.
- **Archivos STEP/STL** (`step_4B3/`):
  - `acoplador_Lac.step` / `.stl`: barra plana 8.86 mm entre agujeros.
  - `poste_Lpc.step` / `.stl`: oreja para Fp, longitud efectiva 8 mm.
  - `poste_Lpd.step` / `.stl`: oreja para Fd, longitud efectiva 18 mm.
  - `pasador.step` / `.stl`: pin de articulacion (D = 1.5 mm).
- **MONTAJE DORSAL** : el 4B#3 se construye del lado DORSAL del dedo (lado
  opuesto a la palma, +Y en el diagrama). Esto es necesario porque las falanges
  apoyan su lado palmar contra los objetos manipulados.
- La reflexion del mecanismo sobre la linea de la bancada (IFP-IFD) preserva
  TODAS las longitudes de eslabon (Lpc=8, Lac=8.86, Lpd=18 mm). La cinematica
  de salida (THETAfd) es identica.

---

## 8. Archivos de referencia

| Archivo | Descripcion |
|---------|-------------|
| `CinematicaExoModificada.m` | Cinematica completa en MATLAB (fuente autoritativa) |
| `diagrama_mecanismo.py` | Genera el diagrama de eslabones general |
| `diagrama_mecanismo_completo.png` | Diagrama completo (postes hsp + 4B#3) |
| `plano_4B3.py` | Genera los planos cotados del 4B#3 |
| `plano_4B3_general.png` | Plano cotado general del 4B#3 montado |
| `plano_4B3_piezas.png` | Plano de fabricacion de las piezas individuales |
| `generar_step_4b3.py` | Genera los STEP/STL de las piezas |
| `step_4B3/*.step` | Archivos CAD de las piezas (formato STEP) |
| `step_4B3/*.stl` | Archivos CAD de las piezas (formato STL) |
| `Analisis_Modificacion_DIP.md` | Documento de analisis del tercer mecanismo |
