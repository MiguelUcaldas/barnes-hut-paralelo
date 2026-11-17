#ifndef ARBOL_H
#define ARBOL_H

#include <stdint.h>

#include "particulas.h"

#define SIN_HIJOS (-1)

/* Las fuentes de una celda ocupan un rango contiguo: primera_particula ..
 * primera_particula + num_particulas - 1.
 *
 * Tres clases de hoja (primer_hijo == SIN_HIJOS):
 *   - vacía:    masa 0, se ignora;
 *   - normal:   num_particulas > 0, se interactúa con cada una;
 *   - resumida: num_particulas == 0 y masa > 0. Solo aparece en árboles
 *               recibidos de otro proceso: es una celda que el dueño no
 *               envió completa porque está lejos de todo este dominio, y
 *               actúa como una sola masa en su centro de masa. */
typedef struct {
    Vector3 centro;            /* centro geométrico del cubo */
    double mitad_lado;
    double masa;
    Vector3 centro_de_masa;
    int primer_hijo;           /* los 8 hijos están seguidos desde aquí */
    int primera_particula;
    int num_particulas;
} Nodo;

/* Arreglo plano de nodos (el 0 es la raíz) y las fuentes de sus hojas. Se
 * usa igual para el árbol local y para los árboles esenciales recibidos. */
typedef struct {
    Nodo *nodos;
    int num_nodos;
    int capacidad;
    const Fuente *fuentes;
} Arbol;

Arbol crear_arbol(int capacidad_inicial);
void liberar_arbol(Arbol *arbol);

/* Construye el árbol de las fuentes locales, ya ordenadas por código Morton.
 * El cubo raíz es el global: es el mismo en todos los procesos, así que las
 * celdas de los árboles de distintos procesos coinciden. */
void construir_arbol_local(Arbol *arbol, const Fuente *fuentes, const uint64_t *codigos, int cantidad,
                           Vector3 centro, double mitad_lado);

#endif
