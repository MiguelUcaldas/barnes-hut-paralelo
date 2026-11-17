#include "particulas.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "memoria.h"

static FILE *abrir_archivo(const char *ruta, const char *modo)
{
    FILE *archivo = fopen(ruta, modo);
    if (archivo == NULL) {
        fprintf(stderr, "Error: no se pudo abrir %s\n", ruta);
        exit(1);
    }
    return archivo;
}

static bool leer_fila(FILE *archivo, double valores[7])
{
    return fscanf(archivo, "%lf %lf %lf %lf %lf %lf %lf",
                  &valores[0], &valores[1], &valores[2], &valores[3],
                  &valores[4], &valores[5], &valores[6]) == 7;
}

/* Primera pasada: contar filas para reservar la memoria exacta de una vez. */
static int contar_filas(FILE *archivo, const char *ruta)
{
    double valores[7];
    int filas = 0;
    while (leer_fila(archivo, valores)) {
        filas++;
    }
    if (!feof(archivo)) {
        fprintf(stderr, "Error: %s, la fila %d no tiene 7 números\n", ruta, filas + 1);
        exit(1);
    }
    rewind(archivo);
    return filas;
}

Particulas crear_particulas(int cantidad)
{
    Particulas particulas;
    particulas.cantidad = cantidad;
    particulas.posiciones = reservar_memoria(cantidad, sizeof(Vector3));
    particulas.velocidades = reservar_memoria(cantidad, sizeof(Vector3));
    particulas.aceleraciones = reservar_memoria(cantidad, sizeof(Vector3));
    particulas.potenciales = reservar_memoria(cantidad, sizeof(double));
    particulas.masas = reservar_memoria(cantidad, sizeof(double));
    particulas.interacciones = reservar_memoria(cantidad, sizeof(int));
    particulas.identificadores = reservar_memoria(cantidad, sizeof(int));
    return particulas;
}

Particulas leer_particulas(const char *ruta)
{
    FILE *archivo = abrir_archivo(ruta, "r");
    int cantidad = contar_filas(archivo, ruta);
    if (cantidad == 0) {
        fprintf(stderr, "Error: %s no tiene partículas\n", ruta);
        exit(1);
    }

    Particulas particulas = crear_particulas(cantidad);
    double valores[7];
    for (int i = 0; i < cantidad; i++) {
        leer_fila(archivo, valores);
        particulas.posiciones[i] = (Vector3){valores[0], valores[1], valores[2]};
        particulas.velocidades[i] = (Vector3){valores[3], valores[4], valores[5]};
        particulas.masas[i] = valores[6];
        particulas.identificadores[i] = i;
        /* Antes del primer cálculo no se conoce el costo: se supone igual
         * para todas. */
        particulas.interacciones[i] = 1;
    }
    fclose(archivo);
    return particulas;
}

/* La partícula que queda en la posición i es la que estaba en nuevo_orden[i].
 * Se copia todo a `auxiliar` y luego se intercambian los arreglos, sin copiar
 * de vuelta. Cada i escribe solo su posición: no hay carreras. */
void reordenar_particulas(Particulas *particulas, Particulas *auxiliar, const int *nuevo_orden)
{
    #pragma omp parallel for
    for (int i = 0; i < particulas->cantidad; i++) {
        int origen = nuevo_orden[i];
        auxiliar->posiciones[i] = particulas->posiciones[origen];
        auxiliar->velocidades[i] = particulas->velocidades[origen];
        auxiliar->aceleraciones[i] = particulas->aceleraciones[origen];
        auxiliar->potenciales[i] = particulas->potenciales[origen];
        auxiliar->masas[i] = particulas->masas[origen];
        auxiliar->interacciones[i] = particulas->interacciones[origen];
        auxiliar->identificadores[i] = particulas->identificadores[origen];
    }
    Particulas intercambio = *particulas;
    *particulas = *auxiliar;
    *auxiliar = intercambio;
}

/* posicion_en_arreglo[fila del archivo] = posición actual en los arreglos.
 * Permite escribir los resultados en el orden del archivo de entrada, que
 * es el que espera verificar.py. */
static int *posiciones_por_fila(const Particulas *particulas)
{
    int *posicion_en_arreglo = reservar_memoria(particulas->cantidad, sizeof(int));
    for (int i = 0; i < particulas->cantidad; i++) {
        posicion_en_arreglo[particulas->identificadores[i]] = i;
    }
    return posicion_en_arreglo;
}

/* Mismo formato y precisión que comun/formato.py, para que verificar.py los lea. */
void escribir_particulas(const char *ruta, const Particulas *particulas)
{
    FILE *archivo = abrir_archivo(ruta, "w");
    int *posicion_en_arreglo = posiciones_por_fila(particulas);
    for (int fila = 0; fila < particulas->cantidad; fila++) {
        int i = posicion_en_arreglo[fila];
        Vector3 posicion = particulas->posiciones[i];
        Vector3 velocidad = particulas->velocidades[i];
        fprintf(archivo, "%.10g %.10g %.10g %.10g %.10g %.10g %.10g\n",
                posicion.x, posicion.y, posicion.z,
                velocidad.x, velocidad.y, velocidad.z, particulas->masas[i]);
    }
    free(posicion_en_arreglo);
    fclose(archivo);
}

void escribir_aceleraciones(const char *ruta, const Particulas *particulas)
{
    FILE *archivo = abrir_archivo(ruta, "w");
    int *posicion_en_arreglo = posiciones_por_fila(particulas);
    for (int fila = 0; fila < particulas->cantidad; fila++) {
        Vector3 aceleracion = particulas->aceleraciones[posicion_en_arreglo[fila]];
        fprintf(archivo, "%.10g %.10g %.10g\n", aceleracion.x, aceleracion.y, aceleracion.z);
    }
    free(posicion_en_arreglo);
    fclose(archivo);
}

void liberar_particulas(Particulas *particulas)
{
    free(particulas->posiciones);
    free(particulas->velocidades);
    free(particulas->aceleraciones);
    free(particulas->potenciales);
    free(particulas->masas);
    free(particulas->interacciones);
    free(particulas->identificadores);
}
