"""Simulación de N cuerpos con Barnes-Hut en Python con multiprocessing.

Cada paso es un leapfrog patada-deriva-patada:
    1. v <- v + a·Δt/2
    2. x <- x + v·Δt
    3. ordenar por Morton, (construir el árbol) y recalcular a en los trabajadores
    4. v <- v + a·Δt/2

Uso:
    python python_mp/barnes_hut.py --entrada datos/uniforme_2000.txt --trabajadores 6
"""

import os

# Cada trabajador usa numpy en un solo hilo. Si numpy abriera además sus
# propios hilos (BLAS), habría más hilos que núcleos y el escalado medido no
# sería el de los procesos. Tiene que ir antes de importar numpy.
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(variable, "1")

import argparse
import re
import time
from pathlib import Path

import numpy as np

from arbol import calcular_cubo_raiz, codigos_morton
from medicion import (TiemposPorFase, agregar_fila_csv, imprimir_encabezado_pasos, imprimir_energia,
                      imprimir_entorno, imprimir_paso, imprimir_resumen)
from trabajadores import EquipoDeTrabajadores


def leer_configuracion():
    parser = argparse.ArgumentParser(description="Barnes-Hut con multiprocessing.", add_help=False)
    parser.add_argument("--ayuda", "-h", action="help", help="muestra esta ayuda")
    parser.add_argument("--entrada", required=True, help="partículas: x y z vx vy vz masa")
    parser.add_argument("--pasos", type=int, default=100)
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--theta", type=float, default=0.5)
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument("--trabajadores", type=int, default=os.cpu_count())
    parser.add_argument("--estrategia", choices=("compartida", "local", "copia"), default="compartida",
                        help="cómo llega el árbol a los trabajadores")
    parser.add_argument("--reparto", choices=("estatico", "dinamico", "costo"), default="dinamico")
    # Con grupos grandes las matrices grupo × fuentes no caben en la caché y
    # varios trabajadores saturan la memoria; 16 fue el mejor tamaño medido.
    parser.add_argument("--tamano-grupo", type=int, default=16,
                        help="partículas que recorren el árbol juntas")
    parser.add_argument("--cada-energia", type=int, default=0,
                        help="imprimir la energía cada K pasos (0 = solo inicio y final)")
    parser.add_argument("--guardar-aceleraciones", help="aceleraciones del estado inicial, para verificar.py")
    parser.add_argument("--guardar-final", help="estado final, para verificar.py")
    parser.add_argument("--distribucion", help="nombre para el CSV (por defecto, sale del archivo)")
    parser.add_argument("--csv", default="resultados/mediciones.csv")
    config = parser.parse_args()

    if config.pasos < 0 or config.dt <= 0 or config.theta < 0 or config.epsilon < 0:
        parser.error("pasos >= 0, dt > 0, theta >= 0 y epsilon >= 0")
    if config.trabajadores < 1 or config.tamano_grupo < 1:
        parser.error("--trabajadores y --tamano-grupo deben ser 1 o más")
    if config.estrategia == "copia" and config.reparto == "dinamico":
        parser.error("con --estrategia copia el árbol viaja dentro de cada tarea; con reparto dinámico "
                     "serían cientos de copias por paso. Use --reparto estatico o costo")
    if config.distribucion is None:
        # "datos/grandes/dos_cumulos_20000.txt" -> "dos_cumulos"
        config.distribucion = re.sub(r"_\d+$", "", Path(config.entrada).stem)
    return config


class Simulacion:
    """Estado de las partículas. Posiciones, masas y aceleraciones están en
    memoria compartida (los trabajadores las leen o escriben); velocidades e
    identificadores solo los usa el proceso principal."""

    def __init__(self, ruta, config):
        datos = np.loadtxt(ruta, ndmin=2)
        if datos.shape[1] != 7:
            raise SystemExit(f"Error: {ruta} debe tener 7 columnas y tiene {datos.shape[1]}")
        self.num_particulas = len(datos)
        self.equipo = EquipoDeTrabajadores(self.num_particulas, config)
        self.compartidas = self.equipo.particulas.arreglos
        self.compartidas["posiciones"][:] = datos[:, 0:3]
        self.compartidas["masas"][:] = datos[:, 6]
        self.velocidades = datos[:, 3:6].copy()
        self.identificadores = np.arange(self.num_particulas)   # fila del archivo de entrada

    @property
    def posiciones(self):
        return self.compartidas["posiciones"]

    @property
    def masas(self):
        return self.compartidas["masas"]

    @property
    def aceleraciones(self):
        return self.compartidas["aceleraciones"]

    def ordenar_por_morton(self, centro, lado):
        """Reacomoda las partículas por código Morton: las de cada celda quedan
        contiguas y los grupos consecutivos quedan cerca en el espacio."""
        codigos = codigos_morton(self.posiciones, centro, lado)
        orden = np.argsort(codigos, kind="stable")
        self.compartidas["codigos"][:] = codigos[orden]
        self.posiciones[:] = self.posiciones[orden]
        self.masas[:] = self.masas[orden]
        self.velocidades = self.velocidades[orden]
        self.identificadores = self.identificadores[orden]

    def en_orden_del_archivo(self, arreglo):
        resultado = np.empty_like(arreglo)
        resultado[self.identificadores] = arreglo
        return resultado

    def energia_total(self):
        cinetica = 0.5 * np.sum(self.masas * np.sum(self.velocidades**2, axis=1))
        # Cada par (i, j) aparece en φ_i y en φ_j, por eso el ½.
        potencial = 0.5 * np.sum(self.masas * self.compartidas["potenciales"])
        return cinetica + potencial


def media_patada(sim, dt):
    sim.velocidades += sim.aceleraciones * (dt / 2)


def deriva(sim, dt):
    sim.posiciones[:] += sim.velocidades * dt


def recalcular_aceleraciones(sim, paso, tiempos):
    """Paso 3 del leapfrog, con la "foto congelada" de las posiciones actuales."""
    marca = time.perf_counter()
    centro, lado = calcular_cubo_raiz(sim.posiciones)
    sim.ordenar_por_morton(centro, lado)
    tiempos.ordenamiento += time.perf_counter() - marca

    resultado = sim.equipo.calcular_aceleraciones(paso, lado)
    tiempos.construccion += resultado["construccion"]
    tiempos.fuerzas += resultado["fuerzas"]
    tiempos.desbalance = resultado["desbalance"]
    tiempos.construccion_local = resultado["construccion_local"]
    return resultado["num_nodos"]


def avanzar_un_paso(sim, config, paso):
    tiempos = TiemposPorFase()
    inicio_paso = time.perf_counter()

    marca = time.perf_counter()
    media_patada(sim, config.dt)
    deriva(sim, config.dt)
    tiempos.integracion = time.perf_counter() - marca

    recalcular_aceleraciones(sim, paso, tiempos)

    marca = time.perf_counter()
    media_patada(sim, config.dt)
    tiempos.integracion += time.perf_counter() - marca

    tiempos.total = time.perf_counter() - inicio_paso
    return tiempos


def ejecutar(sim, config, version):
    imprimir_entorno(version, config, sim.num_particulas)

    # La primera media patada necesita las aceleraciones del estado inicial.
    iniciales = TiemposPorFase()
    num_nodos = recalcular_aceleraciones(sim, 0, iniciales)
    print(f"Cálculo inicial de fuerzas: "
          f"{iniciales.ordenamiento + iniciales.construccion + iniciales.fuerzas:.3f} s ({num_nodos} nodos)\n")
    if config.guardar_aceleraciones:
        np.savetxt(config.guardar_aceleraciones, sim.en_orden_del_archivo(sim.aceleraciones), fmt="%.10g")
    energia_inicial = sim.energia_total()

    acumulado = TiemposPorFase()
    imprimir_encabezado_pasos()
    for paso in range(1, config.pasos + 1):
        tiempos = avanzar_un_paso(sim, config, paso)
        acumulado.sumar(tiempos)
        imprimir_paso(paso, tiempos)
        if config.cada_energia > 0 and paso % config.cada_energia == 0:
            imprimir_energia(paso, sim.energia_total(), energia_inicial)

    energia_final = sim.energia_total()
    imprimir_resumen(acumulado, config.pasos, energia_inicial, energia_final, config.estrategia)
    agregar_fila_csv(config, version, sim.num_particulas, acumulado, energia_inicial, energia_final)

    if config.guardar_final:
        estado = np.column_stack((sim.posiciones, sim.velocidades, sim.masas))
        np.savetxt(config.guardar_final, sim.en_orden_del_archivo(estado), fmt="%.10g")


def main():
    config = leer_configuracion()
    version = f"python_{config.estrategia}_{config.reparto}"
    sim = Simulacion(config.entrada, config)
    try:
        ejecutar(sim, config, version)
    finally:
        # Los bloques de memoria compartida sobreviven al proceso si no se
        # borran: hay que cerrarlos aunque la simulación falle.
        sim.equipo.cerrar()


if __name__ == "__main__":
    main()
