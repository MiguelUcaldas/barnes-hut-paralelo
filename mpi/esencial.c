#include "esencial.h"

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "criterio.h"
#include "memoria.h"

Intercambio crear_intercambio(void)
{
    Intercambio intercambio;
    MPI_Comm_rank(MPI_COMM_WORLD, &intercambio.rango);
    MPI_Comm_size(MPI_COMM_WORLD, &intercambio.num_procesos);
    int p = intercambio.num_procesos;
    intercambio.para = reservar_memoria(p, sizeof(ArbolPodado));
    for (int r = 0; r < p; r++) {
        ArbolPodado *podado = &intercambio.para[r];
        podado->capacidad_nodos = 1024;
        podado->nodos = reservar_memoria(podado->capacidad_nodos, sizeof(Nodo));
        podado->origen = reservar_memoria(podado->capacidad_nodos, sizeof(int));
        podado->capacidad_fuentes = 1024;
        podado->fuentes = reservar_memoria(podado->capacidad_fuentes, sizeof(Fuente));
    }
    intercambio.cajas_minimas = reservar_memoria(p, sizeof(Vector3));
    intercambio.cajas_maximas = reservar_memoria(p, sizeof(Vector3));
    intercambio.nodos_recibidos = NULL;
    intercambio.fuentes_recibidas = NULL;
    intercambio.remotos = reservar_memoria(p, sizeof(Arbol));
    return intercambio;
}

void liberar_intercambio(Intercambio *intercambio)
{
    for (int r = 0; r < intercambio->num_procesos; r++) {
        free(intercambio->para[r].nodos);
        free(intercambio->para[r].origen);
        free(intercambio->para[r].fuentes);
    }
    free(intercambio->para);
    free(intercambio->cajas_minimas);
    free(intercambio->cajas_maximas);
    free(intercambio->nodos_recibidos);
    free(intercambio->fuentes_recibidas);
    free(intercambio->remotos);
}

/* Agrega al árbol podado una copia del nodo local `original` y devuelve su índice. */
static int copiar_nodo(ArbolPodado *podado, const Arbol *local, int original)
{
    if (podado->num_nodos == podado->capacidad_nodos) {
        podado->capacidad_nodos *= 2;
        podado->nodos = ampliar_memoria(podado->nodos, podado->capacidad_nodos, sizeof(Nodo));
        podado->origen = ampliar_memoria(podado->origen, podado->capacidad_nodos, sizeof(int));
    }
    int indice = podado->num_nodos++;
    podado->nodos[indice] = local->nodos[original];
    podado->origen[indice] = original;
    return indice;
}

static int copiar_fuentes(ArbolPodado *podado, const Fuente *fuentes, int cantidad)
{
    while (podado->num_fuentes + cantidad > podado->capacidad_fuentes) {
        podado->capacidad_fuentes *= 2;
        podado->fuentes = ampliar_memoria(podado->fuentes, podado->capacidad_fuentes, sizeof(Fuente));
    }
    int primera = podado->num_fuentes;
    memcpy(&podado->fuentes[primera], fuentes, cantidad * sizeof(Fuente));
    podado->num_fuentes += cantidad;
    return primera;
}

/* Recorre el árbol local con el criterio evaluado contra la caja entera del
 * otro proceso. Se recorre por niveles usando el propio arreglo de salida
 * como cola: los nodos se agregan al final y se procesan en orden. Así los 8
 * hijos de cada nodo quedan seguidos, igual que en el árbol original.
 *   - lejos de toda la caja: se envía resumido (solo masa y centro de masa);
 *   - hoja cercana: se envía con sus partículas;
 *   - nodo interno cercano: se envía y se agregan sus 8 hijos. */
static void extraer_arbol_podado(const Arbol *local, Vector3 caja_minima, Vector3 caja_maxima, double theta,
                                 ArbolPodado *podado)
{
    podado->num_nodos = 0;
    podado->num_fuentes = 0;
    copiar_nodo(podado, local, 0);

    for (int k = 0; k < podado->num_nodos; k++) {
        const Nodo *original = &local->nodos[podado->origen[k]];
        /* Se trabaja con el índice k y no con un puntero a podado->nodos[k]:
         * copiar_nodo puede mover el arreglo con realloc. */
        if (original->masa == 0.0) {
            continue;   /* hoja vacía: se envía tal cual y el receptor la ignora */
        }
        if (celda_suficientemente_lejana(original, caja_minima, caja_maxima, theta)) {
            podado->nodos[k].primer_hijo = SIN_HIJOS;
            podado->nodos[k].num_particulas = 0;
            continue;
        }
        if (original->primer_hijo == SIN_HIJOS) {
            podado->nodos[k].primera_particula = copiar_fuentes(
                podado, &local->fuentes[original->primera_particula], original->num_particulas);
            continue;
        }
        int primer_hijo = podado->num_nodos;
        for (int hijo = 0; hijo < 8; hijo++) {
            copiar_nodo(podado, local, original->primer_hijo + hijo);
        }
        podado->nodos[k].primer_hijo = primer_hijo;
    }
}

static void publicar_cajas(Intercambio *intercambio, const ParticulasLocales *particulas)
{
    double caja[6] = {INFINITY, INFINITY, INFINITY, -INFINITY, -INFINITY, -INFINITY};
    for (int i = 0; i < particulas->cantidad; i++) {
        Vector3 posicion = particulas->datos[i].posicion;
        caja[0] = fmin(caja[0], posicion.x);
        caja[1] = fmin(caja[1], posicion.y);
        caja[2] = fmin(caja[2], posicion.z);
        caja[3] = fmax(caja[3], posicion.x);
        caja[4] = fmax(caja[4], posicion.y);
        caja[5] = fmax(caja[5], posicion.z);
    }
    /* Un proceso sin partículas publica una caja vacía (mínimo +∞, máximo
     * -∞): todo está infinitamente lejos de ella y se le envían solo raíces
     * resumidas. */
    double *todas = reservar_memoria(6 * intercambio->num_procesos, sizeof(double));
    MPI_Allgather(caja, 6, MPI_DOUBLE, todas, 6, MPI_DOUBLE, MPI_COMM_WORLD);
    for (int r = 0; r < intercambio->num_procesos; r++) {
        intercambio->cajas_minimas[r] = (Vector3){todas[6 * r], todas[6 * r + 1], todas[6 * r + 2]};
        intercambio->cajas_maximas[r] = (Vector3){todas[6 * r + 3], todas[6 * r + 4], todas[6 * r + 5]};
    }
    free(todas);
}

/* Desplazamientos (en bytes) a partir de las cantidades de cada proceso. */
static int desplazamientos(const int *bytes, int *desplazamiento, int num_procesos)
{
    int total = 0;
    for (int r = 0; r < num_procesos; r++) {
        desplazamiento[r] = total;
        total += bytes[r];
    }
    return total;
}

/* Envía a cada proceso un bloque de datos distinto con un solo Alltoallv,
 * después de avisar con un Alltoall cuántos bytes recibirá cada uno. Junta
 * en un buffer contiguo lo que va a cada destino, en orden de destino. */
static void *enviar_a_todos(int num_procesos, const int *bytes_a_enviar, void *const *bloques,
                            int *bytes_recibidos, int *desplazamiento_recibido)
{
    int *desplazamiento_envio = reservar_memoria(num_procesos, sizeof(int));
    int total_envio = desplazamientos(bytes_a_enviar, desplazamiento_envio, num_procesos);
    char *envio = reservar_memoria(total_envio > 0 ? total_envio : 1, 1);
    for (int r = 0; r < num_procesos; r++) {
        memcpy(envio + desplazamiento_envio[r], bloques[r], bytes_a_enviar[r]);
    }

    MPI_Alltoall(bytes_a_enviar, 1, MPI_INT, bytes_recibidos, 1, MPI_INT, MPI_COMM_WORLD);
    int total_recibido = desplazamientos(bytes_recibidos, desplazamiento_recibido, num_procesos);
    char *recibido = reservar_memoria(total_recibido > 0 ? total_recibido : 1, 1);
    MPI_Alltoallv(envio, bytes_a_enviar, desplazamiento_envio, MPI_BYTE,
                  recibido, bytes_recibidos, desplazamiento_recibido, MPI_BYTE, MPI_COMM_WORLD);

    free(envio);
    free(desplazamiento_envio);
    return recibido;
}

void intercambiar_arboles_esenciales(Intercambio *intercambio, const Arbol *local, const ParticulasLocales *particulas,
                                     double theta)
{
    int p = intercambio->num_procesos;
    publicar_cajas(intercambio, particulas);

    int *bytes_nodos = reservar_memoria(p, sizeof(int));
    int *bytes_fuentes = reservar_memoria(p, sizeof(int));
    void **bloques_nodos = reservar_memoria(p, sizeof(void *));
    void **bloques_fuentes = reservar_memoria(p, sizeof(void *));
    intercambio->nodos_enviados = 0;
    for (int r = 0; r < p; r++) {
        ArbolPodado *podado = &intercambio->para[r];
        podado->num_nodos = 0;
        podado->num_fuentes = 0;
        if (r != intercambio->rango) {
            extraer_arbol_podado(local, intercambio->cajas_minimas[r], intercambio->cajas_maximas[r], theta, podado);
        }
        bytes_nodos[r] = podado->num_nodos * (int)sizeof(Nodo);
        bytes_fuentes[r] = podado->num_fuentes * (int)sizeof(Fuente);
        bloques_nodos[r] = podado->nodos;
        bloques_fuentes[r] = podado->fuentes;
        intercambio->nodos_enviados += podado->num_nodos;
    }

    int *bytes_nodos_recibidos = reservar_memoria(p, sizeof(int));
    int *desplazamiento_nodos = reservar_memoria(p, sizeof(int));
    int *bytes_fuentes_recibidas = reservar_memoria(p, sizeof(int));
    int *desplazamiento_fuentes = reservar_memoria(p, sizeof(int));
    free(intercambio->nodos_recibidos);
    free(intercambio->fuentes_recibidas);
    intercambio->nodos_recibidos = enviar_a_todos(p, bytes_nodos, bloques_nodos,
                                                  bytes_nodos_recibidos, desplazamiento_nodos);
    intercambio->fuentes_recibidas = enviar_a_todos(p, bytes_fuentes, bloques_fuentes,
                                                    bytes_fuentes_recibidas, desplazamiento_fuentes);

    /* Cada árbol recibido es una vista sobre su parte del buffer. Sus índices
     * (primer_hijo, primera_particula) son relativos a su propio inicio, así
     * que se pueden usar sin traducir. */
    intercambio->nodos_recibidos_total = 0;
    intercambio->fuentes_recibidas_total = 0;
    for (int r = 0; r < p; r++) {
        Arbol *remoto = &intercambio->remotos[r];
        remoto->nodos = (Nodo *)((char *)intercambio->nodos_recibidos + desplazamiento_nodos[r]);
        remoto->num_nodos = bytes_nodos_recibidos[r] / (int)sizeof(Nodo);
        remoto->capacidad = remoto->num_nodos;
        remoto->fuentes = (const Fuente *)((char *)intercambio->fuentes_recibidas + desplazamiento_fuentes[r]);
        intercambio->nodos_recibidos_total += remoto->num_nodos;
        intercambio->fuentes_recibidas_total += bytes_fuentes_recibidas[r] / (int)sizeof(Fuente);
    }

    free(bytes_nodos);
    free(bytes_fuentes);
    free(bloques_nodos);
    free(bloques_fuentes);
    free(bytes_nodos_recibidos);
    free(desplazamiento_nodos);
    free(bytes_fuentes_recibidas);
    free(desplazamiento_fuentes);
}
