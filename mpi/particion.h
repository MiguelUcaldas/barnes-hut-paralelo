#ifndef PARTICION_H
#define PARTICION_H

#include "configuracion.h"
#include "particulas.h"

/* El cubo más pequeño que contiene las partículas de todos los procesos. */
void calcular_cubo_raiz_global(const ParticulasLocales *particulas, Vector3 *centro, double *mitad_lado);

/* Calcula los códigos Morton, decide qué tramo de la curva le toca a cada
 * proceso, envía cada partícula a su dueño y deja las propias ordenadas por
 * código. Es una operación colectiva: todos los procesos deben llamarla. */
void repartir_particulas(ParticulasLocales *particulas, Vector3 centro, double mitad_lado, Particion particion);

#endif
