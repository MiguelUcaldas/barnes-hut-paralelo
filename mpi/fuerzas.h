#ifndef FUERZAS_H
#define FUERZAS_H

#include "arbol.h"
#include "particulas.h"

/* Aceleración, potencial y costo de cada partícula local. La fuente j del
 * árbol local es la partícula j (están en el mismo orden). remotos[r] es el
 * árbol esencial recibido del proceso r; el propio se salta. */
void calcular_aceleraciones(ParticulasLocales *particulas, const Arbol *local, const Arbol *remotos,
                            int num_procesos, int rango, double theta, double epsilon);

/* Sumas de este proceso; hay que sumarlas entre todos para la energía total. */
double energia_cinetica_local(const ParticulasLocales *particulas);
double energia_potencial_local(const ParticulasLocales *particulas);

#endif
