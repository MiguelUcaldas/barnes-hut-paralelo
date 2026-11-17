#ifndef ENTRADA_SALIDA_H
#define ENTRADA_SALIDA_H

#include "particulas.h"

/* El proceso 0 lee el archivo y reparte bloques iguales de filas. En el
 * primer cálculo de fuerzas se vuelven a repartir por la curva de Morton.
 * Devuelve las partículas de este proceso y, en `total`, el N global. */
ParticulasLocales leer_y_repartir(const char *ruta, int *total);

/* Reúnen las partículas en el proceso 0 y escriben en el orden del archivo
 * de entrada, que es el que espera verificar.py. */
void escribir_aceleraciones_globales(const char *ruta, const ParticulasLocales *particulas, int total);
void escribir_estado_global(const char *ruta, const ParticulasLocales *particulas, int total);

#endif
