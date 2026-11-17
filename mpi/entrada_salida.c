#include "entrada_salida.h"

#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "memoria.h"

/* Un error de archivo solo lo ve el proceso 0; MPI_Abort termina a todos,
 * porque los demás estarían esperando en una operación colectiva. */
static void abortar(const char *mensaje, const char *ruta)
{
    fprintf(stderr, "Error: %s %s\n", mensaje, ruta);
    MPI_Abort(MPI_COMM_WORLD, 1);
}

static bool leer_fila(FILE *archivo, double valores[7])
{
    return fscanf(archivo, "%lf %lf %lf %lf %lf %lf %lf",
                  &valores[0], &valores[1], &valores[2], &valores[3],
                  &valores[4], &valores[5], &valores[6]) == 7;
}

static Particula *leer_archivo(const char *ruta, int *total)
{
    FILE *archivo = fopen(ruta, "r");
    if (archivo == NULL) {
        abortar("no se pudo abrir", ruta);
    }
    double valores[7];
    int filas = 0;
    while (leer_fila(archivo, valores)) {
        filas++;
    }
    if (!feof(archivo) || filas == 0) {
        abortar("el archivo no tiene filas de 7 números:", ruta);
    }
    rewind(archivo);

    Particula *particulas = reservar_memoria(filas, sizeof(Particula));
    for (int i = 0; i < filas; i++) {
        leer_fila(archivo, valores);
        particulas[i].posicion = (Vector3){valores[0], valores[1], valores[2]};
        particulas[i].velocidad = (Vector3){valores[3], valores[4], valores[5]};
        particulas[i].masa = valores[6];
        particulas[i].identificador = i;
        /* Antes del primer cálculo no se conoce el costo: igual para todas. */
        particulas[i].costo = 1.0;
    }
    fclose(archivo);
    *total = filas;
    return particulas;
}

ParticulasLocales leer_y_repartir(const char *ruta, int *total)
{
    int rango, num_procesos;
    MPI_Comm_rank(MPI_COMM_WORLD, &rango);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procesos);

    Particula *todas = rango == 0 ? leer_archivo(ruta, total) : NULL;
    MPI_Bcast(total, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int *bytes = reservar_memoria(num_procesos, sizeof(int));
    int *desplazamiento = reservar_memoria(num_procesos, sizeof(int));
    for (int r = 0; r < num_procesos; r++) {
        long desde = (long)*total * r / num_procesos;
        long hasta = (long)*total * (r + 1) / num_procesos;
        bytes[r] = (int)((hasta - desde) * sizeof(Particula));
        desplazamiento[r] = (int)(desde * sizeof(Particula));
    }

    ParticulasLocales locales;
    locales.cantidad = bytes[rango] / (int)sizeof(Particula);
    locales.datos = reservar_memoria(locales.cantidad > 0 ? locales.cantidad : 1, sizeof(Particula));
    MPI_Scatterv(todas, bytes, desplazamiento, MPI_BYTE, locales.datos, bytes[rango], MPI_BYTE, 0, MPI_COMM_WORLD);

    free(todas);
    free(bytes);
    free(desplazamiento);
    return locales;
}

/* Junta todas las partículas en el proceso 0, ordenadas por identificador.
 * En los demás procesos devuelve NULL. */
static Particula *reunir_en_raiz(const ParticulasLocales *particulas, int total)
{
    int rango, num_procesos;
    MPI_Comm_rank(MPI_COMM_WORLD, &rango);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procesos);

    int mis_bytes = particulas->cantidad * (int)sizeof(Particula);
    int *bytes = rango == 0 ? reservar_memoria(num_procesos, sizeof(int)) : NULL;
    MPI_Gather(&mis_bytes, 1, MPI_INT, bytes, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int *desplazamiento = NULL;
    Particula *recibidas = NULL;
    if (rango == 0) {
        desplazamiento = reservar_memoria(num_procesos, sizeof(int));
        for (int r = 1; r < num_procesos; r++) {
            desplazamiento[r] = desplazamiento[r - 1] + bytes[r - 1];
        }
        recibidas = reservar_memoria(total, sizeof(Particula));
    }
    MPI_Gatherv(particulas->datos, mis_bytes, MPI_BYTE, recibidas, bytes, desplazamiento, MPI_BYTE, 0, MPI_COMM_WORLD);

    Particula *en_orden = NULL;
    if (rango == 0) {
        en_orden = reservar_memoria(total, sizeof(Particula));
        for (int i = 0; i < total; i++) {
            en_orden[recibidas[i].identificador] = recibidas[i];
        }
    }
    free(bytes);
    free(desplazamiento);
    free(recibidas);
    return en_orden;
}

void escribir_aceleraciones_globales(const char *ruta, const ParticulasLocales *particulas, int total)
{
    Particula *todas = reunir_en_raiz(particulas, total);
    if (todas == NULL) {
        return;
    }
    FILE *archivo = fopen(ruta, "w");
    if (archivo == NULL) {
        abortar("no se pudo escribir", ruta);
    }
    for (int i = 0; i < total; i++) {
        Vector3 a = todas[i].aceleracion;
        fprintf(archivo, "%.10g %.10g %.10g\n", a.x, a.y, a.z);
    }
    fclose(archivo);
    free(todas);
}

void escribir_estado_global(const char *ruta, const ParticulasLocales *particulas, int total)
{
    Particula *todas = reunir_en_raiz(particulas, total);
    if (todas == NULL) {
        return;
    }
    FILE *archivo = fopen(ruta, "w");
    if (archivo == NULL) {
        abortar("no se pudo escribir", ruta);
    }
    for (int i = 0; i < total; i++) {
        Vector3 x = todas[i].posicion;
        Vector3 v = todas[i].velocidad;
        fprintf(archivo, "%.10g %.10g %.10g %.10g %.10g %.10g %.10g\n", x.x, x.y, x.z, v.x, v.y, v.z, todas[i].masa);
    }
    fclose(archivo);
    free(todas);
}
