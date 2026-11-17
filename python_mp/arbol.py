"""Octree en arreglos planos de numpy, construido nivel por nivel a partir del orden de Morton.

No hay un objeto de Python por nodo: un nodo es un índice k, y cada dato del
nodo está en la posición k de un arreglo (masa[k], centro_de_masa[k], ...).
Así el árbol se construye y se recorre con operaciones de numpy sobre muchos
nodos a la vez, y se puede poner tal cual en memoria compartida.
"""

from dataclasses import dataclass

import numpy as np

# 21 bits por eje: 63 bits en total, caben en un entero de 64 bits. El árbol
# tiene a lo sumo 21 niveles por debajo de la raíz.
BITS_POR_EJE = 21


@dataclass
class Arbol:
    """El nodo 0 es la raíz. Solo existen los nodos con partículas.

    Hijos del nodo k: primer_hijo[k] .. primer_hijo[k] + num_hijos[k] - 1
    (num_hijos[k] = 0 si es hoja). Partículas del nodo k: primera_particula[k]
    .. primera_particula[k] + num_particulas[k] - 1, porque están ordenadas
    por Morton.
    """
    tamano: np.ndarray            # lado del cubo de la celda
    masa: np.ndarray
    centro_de_masa: np.ndarray    # forma (num_nodos, 3)
    primer_hijo: np.ndarray
    num_hijos: np.ndarray
    primera_particula: np.ndarray
    num_particulas: np.ndarray

    @property
    def num_nodos(self):
        return len(self.masa)


def calcular_cubo_raiz(posiciones):
    """El cubo más pequeño que contiene todas las partículas: (centro, lado)."""
    minimo = posiciones.min(axis=0)
    maximo = posiciones.max(axis=0)
    return (minimo + maximo) / 2, float((maximo - minimo).max())


def codigos_morton(posiciones, centro, lado):
    """Divide el cubo en 2^21 rebanadas por eje e intercala los bits de x, y, z.

    Resultado: ... z1 y1 x1 z0 y0 x0. Cada grupo de tres bits es un octante,
    empezando por la raíz en los bits más altos.
    """
    rebanadas_por_eje = 1 << BITS_POR_EJE
    esquina = centro - lado / 2
    # Si todas las partículas están en el mismo punto, el cubo mide 0.
    escala = rebanadas_por_eje / lado if lado > 0 else 0.0
    rebanadas = ((posiciones - esquina) * escala).astype(np.int64)
    # La partícula justo en el borde superior caería en la rebanada 2^21.
    rebanadas = np.clip(rebanadas, 0, rebanadas_por_eje - 1)

    codigos = np.zeros(len(posiciones), dtype=np.int64)
    for bit in range(BITS_POR_EJE):
        for eje in range(3):
            codigos |= ((rebanadas[:, eje] >> bit) & 1) << (3 * bit + eje)
    return codigos


def construir_arbol(codigos, posiciones, masas, lado):
    """Construye el árbol de partículas ya ordenadas por código Morton.

    En el nivel L, las partículas cuyos códigos comparten los primeros 3·L
    bits están en la misma celda, y como están ordenadas, son un tramo
    contiguo del arreglo. Cada nivel se obtiene entonces buscando dónde
    cambia el prefijo, para todas las partículas a la vez. Una celda del
    nivel L es un nodo solo si su padre tenía 2 o más partículas (si tenía
    una, el padre ya era una hoja).
    """
    num_particulas = len(codigos)

    # Por cada nivel: dónde empiezan las partículas de cada nodo nuevo,
    # cuántas tiene y quién es su padre.
    inicios = [np.array([0])]
    cantidades = [np.array([num_particulas])]
    padres = [np.array([-1])]
    niveles = [np.array([0])]
    total_nodos = 1

    # Estado de cada partícula: en qué nodo del nivel actual está y si ese
    # nodo se va a dividir en el nivel siguiente.
    nodo_de_cada_particula = np.zeros(num_particulas, dtype=np.int64)
    se_divide = np.full(num_particulas, num_particulas >= 2)

    for nivel in range(1, BITS_POR_EJE + 1):
        if not se_divide.any():
            break
        prefijos = codigos >> (3 * (BITS_POR_EJE - nivel))
        empieza_celda = np.ones(num_particulas, dtype=bool)
        empieza_celda[1:] = prefijos[1:] != prefijos[:-1]
        inicios_celdas = np.flatnonzero(empieza_celda)
        cantidades_celdas = np.diff(np.append(inicios_celdas, num_particulas))

        es_nodo = se_divide[inicios_celdas]
        nuevos = int(es_nodo.sum())
        ids_celdas = np.full(len(inicios_celdas), -1, dtype=np.int64)
        ids_celdas[es_nodo] = total_nodos + np.arange(nuevos)

        inicios.append(inicios_celdas[es_nodo])
        cantidades.append(cantidades_celdas[es_nodo])
        padres.append(nodo_de_cada_particula[inicios_celdas[es_nodo]])
        niveles.append(np.full(nuevos, nivel))
        total_nodos += nuevos

        nodo_de_cada_particula = np.repeat(ids_celdas, cantidades_celdas)
        se_divide = np.repeat(es_nodo & (cantidades_celdas >= 2), cantidades_celdas)

    inicio = np.concatenate(inicios)
    cantidad = np.concatenate(cantidades)
    padre = np.concatenate(padres)
    nivel = np.concatenate(niveles)

    # Los nodos quedaron ordenados por nivel y, dentro de cada nivel, por
    # padre: los hijos de un mismo padre están seguidos.
    primer_hijo = np.full(total_nodos, -1, dtype=np.int64)
    num_hijos = np.zeros(total_nodos, dtype=np.int64)
    padres_con_hijos, primera_posicion, conteo = np.unique(padre[1:], return_index=True, return_counts=True)
    primer_hijo[padres_con_hijos] = primera_posicion + 1
    num_hijos[padres_con_hijos] = conteo

    # Masa y momento (Σ m·x) de un tramo de partículas como diferencia de
    # sumas acumuladas: una resta por nodo en lugar de un ciclo.
    masa_acumulada = np.concatenate(([0.0], np.cumsum(masas)))
    momento_acumulado = np.vstack((np.zeros((1, 3)), np.cumsum(masas[:, None] * posiciones, axis=0)))
    fin = inicio + cantidad
    masa = masa_acumulada[fin] - masa_acumulada[inicio]
    masa_segura = np.where(masa > 0, masa, 1.0)
    centro_de_masa = (momento_acumulado[fin] - momento_acumulado[inicio]) / masa_segura[:, None]

    return Arbol(
        tamano=lado / 2.0**nivel,
        masa=masa,
        centro_de_masa=centro_de_masa,
        primer_hijo=primer_hijo,
        num_hijos=num_hijos,
        primera_particula=inicio,
        num_particulas=cantidad,
    )
