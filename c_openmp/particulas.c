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

Particulas leer_particulas(const char *ruta)
{
    FILE *archivo = abrir_archivo(ruta, "r");
    Particulas particulas;
    particulas.cantidad = contar_filas(archivo, ruta);
    if (particulas.cantidad == 0) {
        fprintf(stderr, "Error: %s no tiene partículas\n", ruta);
        exit(1);
    }

    int n = particulas.cantidad;
    particulas.posiciones = reservar_memoria(n, sizeof(Vector3));
    particulas.velocidades = reservar_memoria(n, sizeof(Vector3));
    particulas.aceleraciones = reservar_memoria(n, sizeof(Vector3));
    particulas.potenciales = reservar_memoria(n, sizeof(double));
    particulas.masas = reservar_memoria(n, sizeof(double));

    double valores[7];
    for (int i = 0; i < n; i++) {
        leer_fila(archivo, valores);
        particulas.posiciones[i] = (Vector3){valores[0], valores[1], valores[2]};
        particulas.velocidades[i] = (Vector3){valores[3], valores[4], valores[5]};
        particulas.masas[i] = valores[6];
    }
    fclose(archivo);
    return particulas;
}

/* Mismo formato y precisión que comun/formato.py, para que verificar.py los lea. */
void escribir_particulas(const char *ruta, const Particulas *particulas)
{
    FILE *archivo = abrir_archivo(ruta, "w");
    for (int i = 0; i < particulas->cantidad; i++) {
        Vector3 posicion = particulas->posiciones[i];
        Vector3 velocidad = particulas->velocidades[i];
        fprintf(archivo, "%.10g %.10g %.10g %.10g %.10g %.10g %.10g\n",
                posicion.x, posicion.y, posicion.z,
                velocidad.x, velocidad.y, velocidad.z, particulas->masas[i]);
    }
    fclose(archivo);
}

void escribir_aceleraciones(const char *ruta, const Particulas *particulas)
{
    FILE *archivo = abrir_archivo(ruta, "w");
    for (int i = 0; i < particulas->cantidad; i++) {
        Vector3 aceleracion = particulas->aceleraciones[i];
        fprintf(archivo, "%.10g %.10g %.10g\n", aceleracion.x, aceleracion.y, aceleracion.z);
    }
    fclose(archivo);
}

void liberar_particulas(Particulas *particulas)
{
    free(particulas->posiciones);
    free(particulas->velocidades);
    free(particulas->aceleraciones);
    free(particulas->potenciales);
    free(particulas->masas);
}
