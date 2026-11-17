"""Arreglos de numpy en memoria compartida entre procesos.

Cada arreglo vive en un bloque de memoria compartida del sistema operativo.
El proceso que los crea obtiene una descripción (nombre del bloque, forma y
tipo de cada arreglo), que es pequeña y se puede enviar en un mensaje. Los
demás procesos se conectan con esa descripción y ven la misma memoria: nada
se copia.
"""

from multiprocessing import shared_memory

import numpy as np


class ArreglosCompartidos:
    def __init__(self, bloques, arreglos, descripcion, es_dueno):
        self.bloques = bloques
        self.arreglos = arreglos          # nombre -> np.ndarray sobre la memoria compartida
        self.descripcion = descripcion    # nombre -> (nombre del bloque, forma, tipo)
        self.es_dueno = es_dueno

    @classmethod
    def crear(cls, formas):
        """formas: nombre -> (forma, tipo). Los arreglos empiezan en cero."""
        bloques, arreglos, descripcion = [], {}, {}
        for nombre, (forma, tipo) in formas.items():
            tipo = np.dtype(tipo)
            tamano = max(1, int(np.prod(forma)) * tipo.itemsize)
            bloque = shared_memory.SharedMemory(create=True, size=tamano)
            arreglo = np.ndarray(forma, dtype=tipo, buffer=bloque.buf)
            arreglo[...] = 0
            bloques.append(bloque)
            arreglos[nombre] = arreglo
            descripcion[nombre] = (bloque.name, forma, tipo.str)
        return cls(bloques, arreglos, descripcion, es_dueno=True)

    @classmethod
    def conectar(cls, descripcion):
        bloques, arreglos = [], {}
        for nombre, (nombre_bloque, forma, tipo) in descripcion.items():
            # track=False: el dueño es quien borra el bloque al final; si los
            # trabajadores también lo registraran, se intentaría borrar dos veces.
            bloque = shared_memory.SharedMemory(name=nombre_bloque, track=False)
            bloques.append(bloque)
            arreglos[nombre] = np.ndarray(forma, dtype=np.dtype(tipo), buffer=bloque.buf)
        return cls(bloques, arreglos, descripcion, es_dueno=False)

    def cerrar(self):
        # Hay que soltar los arreglos antes de cerrar: mientras un arreglo de
        # numpy apunte a la memoria, el bloque no se puede cerrar.
        self.arreglos.clear()
        for bloque in self.bloques:
            bloque.close()
            if self.es_dueno:
                bloque.unlink()
        self.bloques.clear()
