#ifndef MEMORIA_H
#define MEMORIA_H

#include <stdio.h>
#include <stdlib.h>

/* Si falta memoria no hay forma razonable de seguir: se avisa y se termina,
 * en vez de revisar NULL en cada llamada. */

static inline void *reservar_memoria(size_t cantidad, size_t tamano_elemento)
{
    void *memoria = calloc(cantidad, tamano_elemento);
    if (memoria == NULL) {
        fprintf(stderr, "Error: no hay memoria para %zu elementos de %zu bytes\n", cantidad, tamano_elemento);
        exit(1);
    }
    return memoria;
}

static inline void *ampliar_memoria(void *memoria, size_t cantidad, size_t tamano_elemento)
{
    void *nueva = realloc(memoria, cantidad * tamano_elemento);
    if (nueva == NULL) {
        fprintf(stderr, "Error: no hay memoria para %zu elementos de %zu bytes\n", cantidad, tamano_elemento);
        exit(1);
    }
    return nueva;
}

#endif
