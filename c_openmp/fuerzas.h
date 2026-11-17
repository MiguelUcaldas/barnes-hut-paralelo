#ifndef FUERZAS_H
#define FUERZAS_H

#include "arbol.h"
#include "configuracion.h"
#include "particulas.h"

/* Qué partículas calcula cada hilo y cuánto tardó cada uno. */
typedef struct {
    Reparto reparto;
    int num_hilos;
    int *limites;       /* el hilo h calcula las partículas limites[h] .. limites[h+1]-1 (estático y costo) */
    double *tiempos;    /* segundos de fuerzas de cada hilo en el último cálculo */
} RepartoDeTrabajo;

RepartoDeTrabajo crear_reparto(Reparto reparto, int num_hilos);
void liberar_reparto(RepartoDeTrabajo *reparto);

/* Llena aceleraciones, potenciales e interacciones de todas las partículas.
 * Devuelve el desbalance: tiempo del hilo más lento / tiempo promedio. */
double calcular_aceleraciones(const Arbol *arbol, Particulas *particulas, double theta, double epsilon,
                              RepartoDeTrabajo *reparto);

double energia_cinetica(const Particulas *particulas);
double energia_potencial(const Particulas *particulas);
double energia_total(const Particulas *particulas);

#endif
