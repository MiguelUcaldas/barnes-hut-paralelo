#ifndef ARBOL_H
#define ARBOL_H

#include "particulas.h"

#define SIN_HIJOS (-1)
#define SIN_PARTICULA (-1)

/* A esta profundidad el cubo mide 2^-60 del cubo raíz: dos partículas que
 * siguen juntas ahí están, en la práctica, en el mismo punto. */
#define PROFUNDIDAD_MAXIMA 60

typedef struct {
    Vector3 centro;            /* centro geométrico del cubo */
    double mitad_lado;
    double masa;
    Vector3 centro_de_masa;
    int primer_hijo;           /* los 8 hijos están seguidos desde aquí; SIN_HIJOS si es hoja */
    int particula;             /* en una hoja: su partícula, o SIN_PARTICULA si está vacía */
} Nodo;

/* El árbol es un arreglo plano de nodos; los hijos se referencian por índice.
 * El nodo 0 es la raíz. */
typedef struct {
    Nodo *nodos;
    int num_nodos;
    int capacidad;
    int *siguiente;            /* siguiente[i]: otra partícula en la misma hoja, o SIN_PARTICULA */
    int profundidad_maxima;    /* la alcanzada en la última construcción */
} Arbol;

Arbol crear_arbol(int num_particulas);
void construir_arbol(Arbol *arbol, const Particulas *particulas);
void liberar_arbol(Arbol *arbol);

#endif
