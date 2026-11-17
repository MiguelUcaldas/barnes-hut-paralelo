#ifndef ARBOL_H
#define ARBOL_H

#include <stdint.h>

#include "particulas.h"

#define SIN_HIJOS (-1)

/* Las partículas están ordenadas por Morton, así que las de cualquier celda
 * ocupan un rango contiguo del arreglo: primera_particula ..
 * primera_particula + num_particulas - 1. Una hoja tiene 0 o 1 partículas
 * (más solo si tienen el mismo código Morton, es decir, si están a menos de
 * 2^-21 del lado del cubo raíz). */
typedef struct {
    Vector3 centro;            /* centro geométrico del cubo */
    double mitad_lado;
    double masa;
    Vector3 centro_de_masa;
    int primer_hijo;           /* los 8 hijos están seguidos desde aquí; SIN_HIJOS si es hoja */
    int primera_particula;
    int num_particulas;
} Nodo;

/* Cada hilo reserva nodos de a bloques para no competir por el contador
 * global (ver reservar_hijos). Múltiplo de 8: los hijos van de a 8. */
#define NODOS_POR_BLOQUE 4096

/* Los nodos que el hilo puede usar sin sincronizar: siguiente .. fin-1.
 * El relleno lleva el tamaño a 64 bytes, una línea de caché, para que los
 * bloques de dos hilos no compartan línea (false sharing). */
typedef struct {
    int siguiente;
    int fin;
    char relleno[56];
} BloqueDeNodos;

/* Arreglo plano de nodos; el nodo 0 es la raíz. */
typedef struct {
    Nodo *nodos;
    int num_nodos;           /* nodos reservados, incluido lo que sobró en los bloques */
    int nodos_en_uso;
    int capacidad;
    BloqueDeNodos *bloques;  /* uno por hilo */
    int num_bloques;
} Arbol;

Arbol crear_arbol(int num_particulas);
void construir_arbol(Arbol *arbol, const Particulas *particulas, const uint64_t *codigos,
                     Vector3 centro, double mitad_lado);
void liberar_arbol(Arbol *arbol);

#endif
