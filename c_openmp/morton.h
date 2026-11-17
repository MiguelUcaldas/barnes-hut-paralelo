#ifndef MORTON_H
#define MORTON_H

#include <stdint.h>

#include "particulas.h"

/* 21 bits por eje: 63 bits en total, caben en un entero de 64. El cubo raíz
 * se divide en 2^21 rebanadas por eje, así que el árbol tiene a lo sumo 21
 * niveles por debajo de la raíz. */
#define BITS_POR_EJE 21

typedef struct {
    uint64_t codigo;
    int indice;
} ParticulaConCodigo;

/* Memoria de trabajo del ordenamiento; se reserva una vez y se reutiliza. */
typedef struct {
    ParticulaConCodigo *pares;
    ParticulaConCodigo *auxiliar_pares;
    int *nuevo_orden;
    uint64_t *codigos;        /* códigos ya ordenados: codigos[i] es el de la partícula i */
    Particulas auxiliar;      /* destino de reordenar_particulas */
} OrdenMorton;

OrdenMorton crear_orden_morton(int num_particulas);
void liberar_orden_morton(OrdenMorton *orden);

void calcular_cubo_raiz(const Particulas *particulas, Vector3 *centro, double *mitad_lado);
void ordenar_por_morton(Particulas *particulas, Vector3 centro, double mitad_lado, OrdenMorton *orden);

#endif
