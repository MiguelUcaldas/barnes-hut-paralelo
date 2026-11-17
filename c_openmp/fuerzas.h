#ifndef FUERZAS_H
#define FUERZAS_H

#include "arbol.h"
#include "particulas.h"

/* Llena aceleraciones y potenciales de todas las partículas recorriendo el árbol. */
void calcular_aceleraciones(const Arbol *arbol, Particulas *particulas, double theta, double epsilon);

double energia_cinetica(const Particulas *particulas);
double energia_potencial(const Particulas *particulas);
double energia_total(const Particulas *particulas);

#endif
