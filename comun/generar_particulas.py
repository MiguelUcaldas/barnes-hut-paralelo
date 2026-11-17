"""Genera archivos de partículas para probar las cuatro versiones.

Distribuciones:
    uniforme     posiciones al azar en el cubo unitario y en reposo
                 (la que pide el enunciado para desarrollar).
    dos_cumulos  dos cúmulos densos que se acercan para chocar.
    concentrada  el 90 % de las partículas en un cúmulo pequeño fuera del
                 centro y el resto repartido por todo el cubo.

Las dos últimas castigan un reparto ingenuo del trabajo: si el espacio se
divide en partes iguales, unas pocas partes se quedan con casi todas las
partículas y con casi todo el cálculo.

En todas, cada partícula tiene masa 1/N (masa total 1) y las filas se
barajan, para que ninguna versión dependa del orden del archivo.

Uso:
    python comun/generar_particulas.py --distribucion uniforme --n 2000
    python comun/generar_particulas.py --distribucion dos_cumulos --n 20000 --salida datos/grandes/dos_cumulos_20000.txt
"""

import argparse
from pathlib import Path

import numpy as np

from formato import escribir_particulas

# dos_cumulos
RADIO_CUMULO = 0.04            # radio de escala de Plummer de cada cúmulo
SEPARACION_CUMULOS = 0.4       # distancia inicial entre los dos centros
VELOCIDAD_CHOQUE = 0.5         # rapidez con que se acercan uno al otro

# concentrada
FRACCION_CONCENTRADA = 0.9
RADIO_CONCENTRADO = 0.05
CENTRO_CONCENTRADO = np.array([0.3, 0.3, 0.3])

CENTRO_DEL_CUBO = np.array([0.5, 0.5, 0.5])


def direcciones_al_azar(rng, cantidad):
    """Vectores unitarios repartidos de manera uniforme sobre la esfera.

    Un vector con tres componentes normales no tiene dirección preferida,
    así que al normalizarlo queda un punto uniforme en la esfera.
    """
    vectores = rng.normal(size=(cantidad, 3))
    return vectores / np.linalg.norm(vectores, axis=1, keepdims=True)


def fracciones_de_escape(rng, cantidad):
    """Rapidez de cada partícula como fracción q de la velocidad de escape.

    En el modelo de Plummer, q se distribuye como q² (1 - q²)^(7/2). Se
    muestrea por rechazo: el máximo de esa curva es menor que 0,1.
    """
    aceptadas = np.empty(0)
    while len(aceptadas) < cantidad:
        q = rng.uniform(0.0, 1.0, cantidad)
        altura = rng.uniform(0.0, 0.1, cantidad)
        aceptadas = np.concatenate([aceptadas, q[altura < q**2 * (1.0 - q**2) ** 3.5]])
    return aceptadas[:cantidad]


def esfera_de_plummer(rng, cantidad, masa_cumulo, radio_escala):
    """Cúmulo en equilibrio centrado en el origen (Aarseth, Hénon y Wielen, 1974).

    Las velocidades equilibran la gravedad del propio cúmulo, por eso no
    colapsa de golpe y la conservación de la energía mide el error del
    código y no un colapso violento.
    """
    # Se invierte la masa acumulada M(r)/M = r³ / (r² + a²)^(3/2). Se corta
    # en el 99 % de la masa para no lanzar partículas muy lejos.
    fraccion_de_masa = rng.uniform(0.0, 0.99, cantidad)
    radios = radio_escala / np.sqrt(fraccion_de_masa ** (-2.0 / 3.0) - 1.0)
    posiciones = radios[:, None] * direcciones_al_azar(rng, cantidad)

    velocidad_escape = np.sqrt(2.0 * masa_cumulo / radio_escala) * (1.0 + (radios / radio_escala) ** 2) ** -0.25
    rapideces = fracciones_de_escape(rng, cantidad) * velocidad_escape
    velocidades = rapideces[:, None] * direcciones_al_azar(rng, cantidad)
    return posiciones, velocidades


def uniforme(rng, num_particulas):
    posiciones = rng.random((num_particulas, 3))
    velocidades = np.zeros((num_particulas, 3))
    return posiciones, velocidades


def dos_cumulos(rng, num_particulas):
    cantidad_izquierda = num_particulas // 2
    cantidad_derecha = num_particulas - cantidad_izquierda
    desplazamiento = np.array([SEPARACION_CUMULOS / 2, 0.0, 0.0])
    velocidad_de_cada_uno = np.array([VELOCIDAD_CHOQUE / 2, 0.0, 0.0])

    posiciones_izq, velocidades_izq = esfera_de_plummer(
        rng, cantidad_izquierda, cantidad_izquierda / num_particulas, RADIO_CUMULO)
    posiciones_der, velocidades_der = esfera_de_plummer(
        rng, cantidad_derecha, cantidad_derecha / num_particulas, RADIO_CUMULO)

    posiciones = np.vstack([posiciones_izq + CENTRO_DEL_CUBO - desplazamiento,
                            posiciones_der + CENTRO_DEL_CUBO + desplazamiento])
    velocidades = np.vstack([velocidades_izq + velocidad_de_cada_uno,
                             velocidades_der - velocidad_de_cada_uno])
    return posiciones, velocidades


def concentrada(rng, num_particulas):
    cantidad_cumulo = int(round(FRACCION_CONCENTRADA * num_particulas))
    cantidad_fondo = num_particulas - cantidad_cumulo

    posiciones_cumulo, velocidades_cumulo = esfera_de_plummer(
        rng, cantidad_cumulo, cantidad_cumulo / num_particulas, RADIO_CONCENTRADO)
    posiciones_fondo = rng.random((cantidad_fondo, 3))
    velocidades_fondo = np.zeros((cantidad_fondo, 3))

    posiciones = np.vstack([posiciones_cumulo + CENTRO_CONCENTRADO, posiciones_fondo])
    velocidades = np.vstack([velocidades_cumulo, velocidades_fondo])
    return posiciones, velocidades


DISTRIBUCIONES = {
    "uniforme": uniforme,
    "dos_cumulos": dos_cumulos,
    "concentrada": concentrada,
}


def generar(distribucion, num_particulas, semilla):
    rng = np.random.default_rng(semilla)
    posiciones, velocidades = DISTRIBUCIONES[distribucion](rng, num_particulas)
    orden = rng.permutation(num_particulas)
    masas = np.full(num_particulas, 1.0 / num_particulas)
    return posiciones[orden], velocidades[orden], masas


def main():
    parser = argparse.ArgumentParser(description="Genera un archivo de partículas (x y z vx vy vz masa).")
    parser.add_argument("--distribucion", choices=DISTRIBUCIONES, default="uniforme")
    parser.add_argument("--n", type=int, required=True, help="número de partículas")
    parser.add_argument("--semilla", type=int, default=42)
    parser.add_argument("--salida", help="ruta del archivo (por defecto datos/<distribucion>_<n>.txt)")
    args = parser.parse_args()

    salida = Path(args.salida or f"datos/{args.distribucion}_{args.n}.txt")
    # datos/grandes/ está en .gitignore, así que en un clon nuevo no existe.
    salida.parent.mkdir(parents=True, exist_ok=True)
    posiciones, velocidades, masas = generar(args.distribucion, args.n, args.semilla)
    escribir_particulas(salida, posiciones, velocidades, masas)

    print(f"Distribución: {args.distribucion}, N = {args.n}, semilla = {args.semilla}")
    print(f"Posiciones entre {posiciones.min(axis=0).round(3)} y {posiciones.max(axis=0).round(3)}")
    print(f"Archivo: {salida}")


if __name__ == "__main__":
    main()
