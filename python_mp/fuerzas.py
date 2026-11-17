"""Fuerzas con Barnes-Hut por grupos de partículas, vectorizado con numpy.

Recorrer el árbol partícula por partícula en Python sería muy lento. En vez
de eso, se toma un grupo de partículas consecutivas en orden Morton (que
están cerca en el espacio) y se recorre el árbol una vez para todo el
grupo. El resultado es una lista de interacción: celdas lejanas (cada una
actúa como una partícula en su centro de masa) y partículas de las hojas
cercanas. Después se calcula de una vez la matriz grupo × lista.
"""

import numpy as np


def expandir_rangos(primeros, cantidades):
    """Todos los enteros de los rangos [primeros[k], primeros[k] + cantidades[k]).

    Ejemplo: primeros = [10, 50], cantidades = [3, 2] -> [10, 11, 12, 50, 51].
    """
    total = int(cantidades.sum())
    desplazamientos = np.arange(total) - np.repeat(np.cumsum(cantidades) - cantidades, cantidades)
    return np.repeat(primeros, cantidades) + desplazamientos


# CRITERIO DE APERTURA: aquí se decide qué celdas se tratan como una sola
# partícula en su centro de masa (las que devuelven True) y cuáles se abren.
# Enunciado: tamaño / distancia < θ, comparado con cuadrados para no sacar raíces.
# La distancia se mide desde el centro de masa de la celda hasta el punto más
# cercano de la caja del grupo: si la celda está lejos de ese punto, está lejos
# de todas las partículas del grupo. Con θ = 0 nunca se cumple: suma directa.
def celdas_suficientemente_lejanas(arbol, nodos, caja_minima, caja_maxima, theta):
    centros = arbol.centro_de_masa[nodos]
    separacion = np.maximum(caja_minima - centros, 0.0) + np.maximum(centros - caja_maxima, 0.0)
    distancia2 = np.einsum("ij,ij->i", separacion, separacion)
    return arbol.tamano[nodos] ** 2 < theta**2 * distancia2


def lista_de_interaccion(arbol, caja_minima, caja_maxima, theta):
    """Recorre el árbol nivel por nivel para un grupo.

    frontera tiene todos los nodos de un nivel que hay que revisar. Los
    lejanos van a la lista, las hojas cercanas aportan sus partículas, y los
    nodos internos cercanos se abren: sus hijos son la frontera siguiente.
    Devuelve (índices de celdas lejanas, índices de partículas).
    """
    frontera = np.array([0])
    celdas = []
    hojas = []
    while len(frontera) > 0:
        lejanas = celdas_suficientemente_lejanas(arbol, frontera, caja_minima, caja_maxima, theta)
        es_hoja = arbol.num_hijos[frontera] == 0
        celdas.append(frontera[lejanas])
        hojas.append(frontera[~lejanas & es_hoja])
        abrir = frontera[~lejanas & ~es_hoja]
        frontera = expandir_rangos(arbol.primer_hijo[abrir], arbol.num_hijos[abrir])

    celdas = np.concatenate(celdas)
    hojas = np.concatenate(hojas)
    particulas = expandir_rangos(arbol.primera_particula[hojas], arbol.num_particulas[hojas])
    return celdas, particulas


def fuerzas_de_grupo(arbol, posiciones, masas, primera, ultima, theta, epsilon2):
    """Aceleración y potencial de las partículas primera .. ultima - 1.

    Devuelve (aceleraciones, potenciales, interacciones calculadas).
    """
    grupo = posiciones[primera:ultima]
    celdas, particulas = lista_de_interaccion(arbol, grupo.min(axis=0), grupo.max(axis=0), theta)
    fuentes = np.concatenate((arbol.centro_de_masa[celdas], posiciones[particulas]))
    masas_fuente = np.concatenate((arbol.masa[celdas], masas[particulas]))

    # Matrices (partículas del grupo) × (fuentes); un eje a la vez porque los
    # arreglos 2D contiguos son más rápidos que uno 3D.
    dx = fuentes[:, 0] - grupo[:, 0, None]
    dy = fuentes[:, 1] - grupo[:, 1, None]
    dz = fuentes[:, 2] - grupo[:, 2, None]
    with np.errstate(divide="ignore"):
        inversa = 1.0 / np.sqrt(dx * dx + dy * dy + dz * dz + epsilon2)

    # Las partículas del grupo también están en la lista (sus hojas nunca
    # están lejos del grupo). Una partícula no se atrae a sí misma.
    propias = np.flatnonzero((particulas >= primera) & (particulas < ultima))
    inversa[particulas[propias] - primera, len(celdas) + propias] = 0.0

    peso = masas_fuente * (inversa * inversa * inversa)
    aceleraciones = np.column_stack((
        np.einsum("gk,gk->g", peso, dx),
        np.einsum("gk,gk->g", peso, dy),
        np.einsum("gk,gk->g", peso, dz),
    ))
    potenciales = -(inversa @ masas_fuente)
    interacciones = inversa.size - len(propias)
    return aceleraciones, potenciales, interacciones
