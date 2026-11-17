"""Procesos trabajadores, cómo les llega el árbol y cómo se reparte el trabajo.

En cada paso el proceso principal ordena las partículas, reparte los grupos
entre los trabajadores y espera a que todos terminen (Pool.map no regresa
hasta tener todos los resultados: es la barrera de la "foto congelada").
Cada trabajador escribe las aceleraciones de sus grupos directamente en
memoria compartida; como cada grupo es un tramo distinto del arreglo, dos
trabajadores nunca escriben en el mismo lugar.

Estrategias para que el árbol llegue a los trabajadores (--estrategia):
    compartida  el principal lo construye en memoria compartida y los
                trabajadores lo leen sin copiarlo.
    local       cada trabajador construye su propia copia a partir de las
                partículas compartidas (trabajo repetido, cero comunicación).
    copia       el principal lo construye y lo envía serializado (pickle)
                dentro de cada tarea (comunicación explícita).
En las tres, las partículas y los resultados están en memoria compartida:
lo único que cambia es cómo llega el árbol.
"""

import multiprocessing
import os
import time
from dataclasses import dataclass, fields

import numpy as np

from arbol import Arbol, construir_arbol
from fuerzas import fuerzas_de_grupo
from memoria_compartida import ArreglosCompartidos

# En el reparto dinámico, cada tarea tiene este número de grupos. Más pequeño
# balancea mejor, pero cada tarea cuesta un mensaje de ida y uno de vuelta.
GRUPOS_POR_TAREA = 4


@dataclass
class Tarea:
    paso: int
    primer_grupo: int
    fin_grupo: int
    tamano_grupo: int
    theta: float
    epsilon: float
    estrategia: str
    lado: float                              # para construir el árbol (estrategia local)
    descripcion_arbol: dict | None = None    # estrategia compartida
    num_nodos: int = 0                       # estrategia compartida
    arbol: Arbol | None = None               # estrategia copia


# Estado de cada proceso trabajador. Es una variable global del proceso, no
# compartida: cada trabajador tiene la suya.
_trabajador = {}


def iniciar_trabajador(descripcion_particulas):
    """Se ejecuta una vez al crear cada trabajador: se conecta a las partículas."""
    _trabajador["particulas"] = ArreglosCompartidos.conectar(descripcion_particulas)
    _trabajador["arbol_compartido"] = None
    _trabajador["arbol_local"] = None
    _trabajador["paso_del_arbol_local"] = None


def arbol_desde_memoria_compartida(descripcion, num_nodos):
    actual = _trabajador["arbol_compartido"]
    # El principal crea un bloque nuevo cuando el árbol ya no cabe en el
    # anterior: si cambió el nombre, hay que reconectarse.
    if actual is None or actual.descripcion != descripcion:
        if actual is not None:
            actual.cerrar()
        actual = ArreglosCompartidos.conectar(descripcion)
        _trabajador["arbol_compartido"] = actual
    return Arbol(**{nombre: arreglo[:num_nodos] for nombre, arreglo in actual.arreglos.items()})


def obtener_arbol(tarea, particulas):
    """Devuelve (árbol, segundos que se gastaron construyéndolo en este trabajador)."""
    if tarea.estrategia == "copia":
        return tarea.arbol, 0.0
    if tarea.estrategia == "compartida":
        return arbol_desde_memoria_compartida(tarea.descripcion_arbol, tarea.num_nodos), 0.0
    # local: el primer bloque de cada paso que le toca a este trabajador
    # construye el árbol; los demás bloques del mismo paso lo reutilizan.
    if _trabajador["paso_del_arbol_local"] != tarea.paso:
        inicio = time.perf_counter()
        _trabajador["arbol_local"] = construir_arbol(
            particulas["codigos"], particulas["posiciones"], particulas["masas"], tarea.lado)
        _trabajador["paso_del_arbol_local"] = tarea.paso
        return _trabajador["arbol_local"], time.perf_counter() - inicio
    return _trabajador["arbol_local"], 0.0


def calcular_tarea(tarea):
    """Fuerzas de los grupos primer_grupo .. fin_grupo - 1 de una tarea."""
    inicio = time.perf_counter()
    particulas = _trabajador["particulas"].arreglos
    arbol, segundos_construccion = obtener_arbol(tarea, particulas)

    epsilon2 = tarea.epsilon**2
    num_particulas = len(particulas["masas"])
    for grupo in range(tarea.primer_grupo, tarea.fin_grupo):
        primera = grupo * tarea.tamano_grupo
        ultima = min(primera + tarea.tamano_grupo, num_particulas)
        aceleraciones, potenciales, interacciones = fuerzas_de_grupo(
            arbol, particulas["posiciones"], particulas["masas"], primera, ultima, tarea.theta, epsilon2)
        particulas["aceleraciones"][primera:ultima] = aceleraciones
        particulas["potenciales"][primera:ultima] = potenciales
        particulas["costos_de_grupo"][grupo] = interacciones

    return os.getpid(), time.perf_counter() - inicio, segundos_construccion, arbol.num_nodos


# REPARTO DE TRABAJO ENTRE PROCESOS: aquí se decide qué grupos calcula cada
# trabajador. Devuelve una lista de tareas (primer_grupo, fin_grupo).
#   estatico  un tramo por trabajador, con la misma cantidad de grupos.
#   dinamico  tramos de GRUPOS_POR_TAREA; el Pool le da el siguiente al
#             trabajador que termine primero.
#   costo     un tramo por trabajador, con la misma suma de interacciones
#             del paso anterior.
# Los grupos están en orden Morton: un tramo contiguo es una región compacta.
def repartir_grupos(num_grupos, num_trabajadores, reparto, costos_de_grupo):
    if reparto == "dinamico":
        return [(primero, min(primero + GRUPOS_POR_TAREA, num_grupos))
                for primero in range(0, num_grupos, GRUPOS_POR_TAREA)]
    if reparto == "estatico":
        limites = np.linspace(0, num_grupos, num_trabajadores + 1).round().astype(int)
    else:
        costo_acumulado = np.cumsum(costos_de_grupo)
        objetivos = costo_acumulado[-1] * np.arange(1, num_trabajadores) / num_trabajadores
        cortes = np.minimum(np.searchsorted(costo_acumulado, objetivos) + 1, num_grupos)
        limites = np.concatenate(([0], cortes, [num_grupos]))
    return [(int(limites[k]), int(limites[k + 1]))
            for k in range(num_trabajadores) if limites[k] < limites[k + 1]]


class EquipoDeTrabajadores:
    """Los procesos trabajadores y la memoria compartida con ellos."""

    def __init__(self, num_particulas, config):
        self.config = config
        self.num_grupos = -(-num_particulas // config.tamano_grupo)   # división hacia arriba
        self.particulas = ArreglosCompartidos.crear({
            "posiciones": ((num_particulas, 3), np.float64),
            "masas": ((num_particulas,), np.float64),
            "codigos": ((num_particulas,), np.int64),
            "aceleraciones": ((num_particulas, 3), np.float64),
            "potenciales": ((num_particulas,), np.float64),
            "costos_de_grupo": ((self.num_grupos,), np.float64),
        })
        # Antes del primer cálculo no se conoce el costo: igual para todos.
        self.particulas.arreglos["costos_de_grupo"][:] = 1.0
        self.arbol_compartido = None

        # forkserver: cada trabajador arranca desde un proceso limpio y no
        # hereda la memoria del principal. Todo lo que comparten es explícito.
        contexto = multiprocessing.get_context("forkserver")
        self.pool = contexto.Pool(config.trabajadores, initializer=iniciar_trabajador,
                                  initargs=(self.particulas.descripcion,))

    def publicar_arbol(self, arbol):
        """Copia el árbol a memoria compartida; si no cabe, crea un bloque más grande."""
        if self.arbol_compartido is None or len(self.arbol_compartido.arreglos["masa"]) < arbol.num_nodos:
            if self.arbol_compartido is not None:
                self.arbol_compartido.cerrar()
            capacidad = 2 * arbol.num_nodos
            formas = {}
            for campo in fields(Arbol):
                arreglo = getattr(arbol, campo.name)
                formas[campo.name] = ((capacidad,) + arreglo.shape[1:], arreglo.dtype)
            self.arbol_compartido = ArreglosCompartidos.crear(formas)
        for nombre, destino in self.arbol_compartido.arreglos.items():
            destino[:arbol.num_nodos] = getattr(arbol, nombre)

    def calcular_aceleraciones(self, paso, lado):
        """Llena aceleraciones, potenciales y costos en la memoria compartida.

        Devuelve un diccionario con los tiempos de construcción (en el
        principal), de fuerzas (todo el Pool.map), la construcción más lenta
        dentro de un trabajador, el desbalance y el número de nodos.
        """
        config = self.config
        particulas = self.particulas.arreglos
        resultado = {"construccion": 0.0}

        arbol = None
        if config.estrategia != "local":
            inicio = time.perf_counter()
            arbol = construir_arbol(particulas["codigos"], particulas["posiciones"], particulas["masas"], lado)
            if config.estrategia == "compartida":
                self.publicar_arbol(arbol)
            resultado["construccion"] = time.perf_counter() - inicio

        tareas = []
        for primero, fin in repartir_grupos(self.num_grupos, config.trabajadores, config.reparto,
                                            particulas["costos_de_grupo"]):
            tarea = Tarea(paso, primero, fin, config.tamano_grupo, config.theta, config.epsilon,
                          config.estrategia, lado)
            if config.estrategia == "compartida":
                tarea.descripcion_arbol = self.arbol_compartido.descripcion
                tarea.num_nodos = arbol.num_nodos
            elif config.estrategia == "copia":
                tarea.arbol = arbol
            tareas.append(tarea)

        inicio = time.perf_counter()
        respuestas = self.pool.map(calcular_tarea, tareas, chunksize=1)
        resultado["fuerzas"] = time.perf_counter() - inicio

        tiempo_por_trabajador = {}
        for pid, segundos, _, _ in respuestas:
            tiempo_por_trabajador[pid] = tiempo_por_trabajador.get(pid, 0.0) + segundos
        promedio = sum(tiempo_por_trabajador.values()) / config.trabajadores
        resultado["desbalance"] = max(tiempo_por_trabajador.values()) / promedio if promedio > 0 else 1.0
        resultado["construccion_local"] = max(respuesta[2] for respuesta in respuestas)
        resultado["num_nodos"] = respuestas[0][3]
        return resultado

    def cerrar(self):
        self.pool.close()
        self.pool.join()
        if self.arbol_compartido is not None:
            self.arbol_compartido.cerrar()
        self.particulas.cerrar()
