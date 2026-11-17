#include "particion.h"

#include <math.h>
#include <mpi.h>
#include <stdlib.h>

#include "memoria.h"
#include "morton.h"

/* Cada proceso aporta este número de muestras para decidir los cortes. Con
 * 256 muestras por proceso, el error de balance por muestreo queda por
 * debajo del 1 % del trabajo de un proceso. */
#define MUESTRAS_POR_PROCESO 256

typedef struct {
    uint64_t codigo;
    double peso;
} Muestra;

void calcular_cubo_raiz_global(const ParticulasLocales *particulas, Vector3 *centro, double *mitad_lado)
{
    double minimo[3] = {INFINITY, INFINITY, INFINITY};
    double maximo[3] = {-INFINITY, -INFINITY, -INFINITY};
    for (int i = 0; i < particulas->cantidad; i++) {
        Vector3 posicion = particulas->datos[i].posicion;
        minimo[0] = fmin(minimo[0], posicion.x);
        minimo[1] = fmin(minimo[1], posicion.y);
        minimo[2] = fmin(minimo[2], posicion.z);
        maximo[0] = fmax(maximo[0], posicion.x);
        maximo[1] = fmax(maximo[1], posicion.y);
        maximo[2] = fmax(maximo[2], posicion.z);
    }
    /* Cada proceso conoce solo sus partículas: el mínimo y el máximo
     * globales se combinan entre todos y todos reciben el resultado. */
    MPI_Allreduce(MPI_IN_PLACE, minimo, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, maximo, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    *centro = (Vector3){(minimo[0] + maximo[0]) / 2, (minimo[1] + maximo[1]) / 2, (minimo[2] + maximo[2]) / 2};
    *mitad_lado = fmax(maximo[0] - minimo[0], fmax(maximo[1] - minimo[1], maximo[2] - minimo[2])) / 2;
}

static void calcular_codigos(ParticulasLocales *particulas, Vector3 centro, double mitad_lado)
{
    Vector3 esquina = {centro.x - mitad_lado, centro.y - mitad_lado, centro.z - mitad_lado};
    /* Si todas las partículas están en el mismo punto, el cubo mide 0. */
    double escala = mitad_lado > 0.0 ? (1u << BITS_POR_EJE) / (2.0 * mitad_lado) : 0.0;
    for (int i = 0; i < particulas->cantidad; i++) {
        particulas->datos[i].codigo = codigo_morton(particulas->datos[i].posicion, esquina, escala);
    }
}

/* El identificador desempata: el orden es siempre el mismo. */
static int comparar_particulas(const void *a, const void *b)
{
    const Particula *primera = a;
    const Particula *segunda = b;
    if (primera->codigo != segunda->codigo) {
        return primera->codigo < segunda->codigo ? -1 : 1;
    }
    return (primera->identificador > segunda->identificador) - (primera->identificador < segunda->identificador);
}

static int comparar_muestras(const void *a, const void *b)
{
    const Muestra *primera = a;
    const Muestra *segunda = b;
    return (primera->codigo > segunda->codigo) - (primera->codigo < segunda->codigo);
}

static double peso(const Particula *particula, Particion particion)
{
    return particion == PARTICION_COSTO ? particula->costo : 1.0;
}

/* Muestras repartidas a igual peso a lo largo de las partículas locales (ya
 * ordenadas): cada una representa 1/256 del peso de este proceso. */
static void tomar_muestras(const ParticulasLocales *particulas, Particion particion, Muestra *muestras)
{
    int cantidad = particulas->cantidad;
    double total = 0.0;
    for (int i = 0; i < cantidad; i++) {
        total += peso(&particulas->datos[i], particion);
    }
    double peso_de_muestra = total / MUESTRAS_POR_PROCESO;

    double acumulado = 0.0;   /* peso de las partículas antes de la i */
    int i = 0;
    for (int s = 0; s < MUESTRAS_POR_PROCESO; s++) {
        double objetivo = (s + 0.5) * peso_de_muestra;
        while (i < cantidad - 1 && acumulado + peso(&particulas->datos[i], particion) < objetivo) {
            acumulado += peso(&particulas->datos[i], particion);
            i++;
        }
        muestras[s].codigo = cantidad > 0 ? particulas->datos[i].codigo : 0;
        muestras[s].peso = cantidad > 0 ? peso_de_muestra : 0.0;
    }
}

/* PARTICIONADO: aquí se decide qué tramo de la curva de Morton le toca a
 * cada proceso. El proceso k recibe los códigos desde cortes[k-1] (o desde 0
 * si k = 0) hasta antes de cortes[k] (o hasta el final si k es el último).
 * Un tramo contiguo de la curva es una región compacta del espacio.
 *   volumen:  tramos de igual longitud en la curva, es decir, de igual
 *             volumen; no mira dónde están las partículas (reparto ingenuo).
 *   cantidad: cortes a igual número de partículas.
 *   costo:    cortes a igual número de interacciones del paso anterior.
 * Para cantidad y costo, cada proceso aporta muestras de sus partículas, y
 * todos reciben todas las muestras (Allgather) y calculan los mismos cortes. */
static void calcular_cortes(const ParticulasLocales *particulas, Particion particion, int num_procesos,
                            uint64_t *cortes)
{
    if (particion == PARTICION_VOLUMEN) {
        uint64_t longitud_de_la_curva = (uint64_t)1 << (3 * BITS_POR_EJE);
        for (int k = 1; k < num_procesos; k++) {
            cortes[k - 1] = longitud_de_la_curva / num_procesos * k;
        }
        return;
    }

    Muestra locales[MUESTRAS_POR_PROCESO];
    tomar_muestras(particulas, particion, locales);
    int num_muestras = num_procesos * MUESTRAS_POR_PROCESO;
    Muestra *todas = reservar_memoria(num_muestras, sizeof(Muestra));
    MPI_Allgather(locales, MUESTRAS_POR_PROCESO * sizeof(Muestra), MPI_BYTE,
                  todas, MUESTRAS_POR_PROCESO * sizeof(Muestra), MPI_BYTE, MPI_COMM_WORLD);
    qsort(todas, num_muestras, sizeof(Muestra), comparar_muestras);

    double total = 0.0;
    for (int s = 0; s < num_muestras; s++) {
        total += todas[s].peso;
    }
    double acumulado = 0.0;
    int s = 0;
    for (int k = 1; k < num_procesos; k++) {
        double objetivo = total * k / num_procesos;
        while (s < num_muestras - 1 && acumulado + todas[s].peso < objetivo) {
            acumulado += todas[s].peso;
            s++;
        }
        cortes[k - 1] = todas[s].codigo;
    }
    free(todas);
}

/* Envía a cada proceso las partículas de su tramo. Antes de mover los datos,
 * cada proceso tiene que saber cuántas va a recibir de cada uno (Alltoall de
 * cantidades); después se mueven todas con un solo Alltoallv. */
static void intercambiar_particulas(ParticulasLocales *particulas, const uint64_t *cortes, int num_procesos)
{
    int *bytes_a_enviar = reservar_memoria(num_procesos, sizeof(int));
    int *bytes_a_recibir = reservar_memoria(num_procesos, sizeof(int));
    int *desplazamiento_envio = reservar_memoria(num_procesos, sizeof(int));
    int *desplazamiento_recibo = reservar_memoria(num_procesos, sizeof(int));

    /* Las partículas están ordenadas, así que las de cada destino van seguidas. */
    int destino = 0;
    for (int i = 0; i < particulas->cantidad; i++) {
        while (destino < num_procesos - 1 && particulas->datos[i].codigo >= cortes[destino]) {
            destino++;
        }
        bytes_a_enviar[destino] += sizeof(Particula);
    }
    MPI_Alltoall(bytes_a_enviar, 1, MPI_INT, bytes_a_recibir, 1, MPI_INT, MPI_COMM_WORLD);

    int total_a_recibir = 0;
    for (int r = 0; r < num_procesos; r++) {
        desplazamiento_envio[r] = r == 0 ? 0 : desplazamiento_envio[r - 1] + bytes_a_enviar[r - 1];
        desplazamiento_recibo[r] = r == 0 ? 0 : desplazamiento_recibo[r - 1] + bytes_a_recibir[r - 1];
        total_a_recibir += bytes_a_recibir[r];
    }

    int cantidad_nueva = total_a_recibir / (int)sizeof(Particula);
    Particula *recibidas = reservar_memoria(cantidad_nueva > 0 ? cantidad_nueva : 1, sizeof(Particula));
    MPI_Alltoallv(particulas->datos, bytes_a_enviar, desplazamiento_envio, MPI_BYTE,
                  recibidas, bytes_a_recibir, desplazamiento_recibo, MPI_BYTE, MPI_COMM_WORLD);

    free(particulas->datos);
    particulas->datos = recibidas;
    particulas->cantidad = cantidad_nueva;
    free(bytes_a_enviar);
    free(bytes_a_recibir);
    free(desplazamiento_envio);
    free(desplazamiento_recibo);
}

void repartir_particulas(ParticulasLocales *particulas, Vector3 centro, double mitad_lado, Particion particion)
{
    int num_procesos;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procesos);

    calcular_codigos(particulas, centro, mitad_lado);
    qsort(particulas->datos, particulas->cantidad, sizeof(Particula), comparar_particulas);

    uint64_t *cortes = reservar_memoria(num_procesos, sizeof(uint64_t));
    calcular_cortes(particulas, particion, num_procesos, cortes);
    intercambiar_particulas(particulas, cortes, num_procesos);
    free(cortes);

    /* Llegan tramos ordenados de varios procesos, intercalados entre sí. */
    qsort(particulas->datos, particulas->cantidad, sizeof(Particula), comparar_particulas);
}
