"""Lectura y escritura de los archivos que comparten las cuatro versiones.

Partículas: una fila por partícula, siete columnas separadas por espacios:
    x y z vx vy vz masa
Aceleraciones: una fila por partícula, tres columnas, en el mismo orden
que el archivo de partículas del que salieron:
    ax ay az
"""

import numpy as np

# Diez cifras significativas bastan: el error de redondeo (1e-10) queda muy
# por debajo de los errores que se miden (Barnes-Hut, integrador).
FORMATO_NUMERO = "%.10g"


def leer_particulas(ruta):
    datos = np.loadtxt(ruta, ndmin=2)
    if datos.shape[1] != 7:
        raise ValueError(f"{ruta}: se esperaban 7 columnas y hay {datos.shape[1]}")
    posiciones = datos[:, 0:3].copy()
    velocidades = datos[:, 3:6].copy()
    masas = datos[:, 6].copy()
    return posiciones, velocidades, masas


def escribir_particulas(ruta, posiciones, velocidades, masas):
    datos = np.column_stack([posiciones, velocidades, masas])
    np.savetxt(ruta, datos, fmt=FORMATO_NUMERO)


def leer_aceleraciones(ruta):
    datos = np.loadtxt(ruta, ndmin=2)
    if datos.shape[1] != 3:
        raise ValueError(f"{ruta}: se esperaban 3 columnas y hay {datos.shape[1]}")
    return datos


def escribir_aceleraciones(ruta, aceleraciones):
    np.savetxt(ruta, aceleraciones, fmt=FORMATO_NUMERO)
