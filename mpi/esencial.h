#ifndef ESENCIAL_H
#define ESENCIAL_H

#include "arbol.h"
#include "particulas.h"

/* Un árbol podado: la parte del árbol local que otro proceso necesita. */
typedef struct {
    Nodo *nodos;
    int num_nodos;
    int capacidad_nodos;
    Fuente *fuentes;
    int num_fuentes;
    int capacidad_fuentes;
    int *origen;               /* origen[k]: índice en el árbol local del nodo k */
} ArbolPodado;

/* Todo lo que se reutiliza de un paso al siguiente en el intercambio. */
typedef struct {
    int rango;
    int num_procesos;
    ArbolPodado *para;         /* para[r]: lo que se le envía al proceso r */
    Vector3 *cajas_minimas;    /* caja del dominio de cada proceso */
    Vector3 *cajas_maximas;
    Nodo *nodos_recibidos;
    Fuente *fuentes_recibidas;
    Arbol *remotos;            /* remotos[r]: el árbol esencial recibido del proceso r */
    long nodos_enviados;       /* totales del último intercambio, para el resumen */
    long nodos_recibidos_total;
    long fuentes_recibidas_total;
} Intercambio;

Intercambio crear_intercambio(void);
void liberar_intercambio(Intercambio *intercambio);

/* Cada proceso publica la caja de su dominio, extrae para cada otro proceso
 * la parte de su árbol que ese proceso necesita, y la envía. Al terminar,
 * intercambio->remotos[r] es el árbol esencial recibido del proceso r.
 * Operación colectiva: todos los procesos deben llamarla. */
void intercambiar_arboles_esenciales(Intercambio *intercambio, const Arbol *local, const ParticulasLocales *particulas,
                                     double theta);

#endif
