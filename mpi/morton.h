#ifndef MORTON_H
#define MORTON_H

#include <stdint.h>

#include "particulas.h"

/* 21 bits por eje: 63 bits en total. El árbol tiene a lo sumo 21 niveles. */
#define BITS_POR_EJE 21

/* Código Morton de una posición dentro del cubo raíz (esquina inferior y
 * rebanadas por unidad de longitud). */
uint64_t codigo_morton(Vector3 posicion, Vector3 esquina, double escala);

/* Los tres bits del código que dicen en qué hijo cae, en un nodo de este nivel. */
int octante_en_nivel(uint64_t codigo, int nivel);

#endif
