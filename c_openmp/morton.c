#include "morton.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "memoria.h"

/* Por debajo de este tamaño, un bloque se ordena con qsort en un solo hilo:
 * crear más tareas costaría más que lo que se gana. */
#define UMBRAL_ORDENAMIENTO_SECUENCIAL 8192

OrdenMorton crear_orden_morton(int num_particulas)
{
    OrdenMorton orden;
    orden.pares = reservar_memoria(num_particulas, sizeof(ParticulaConCodigo));
    orden.auxiliar_pares = reservar_memoria(num_particulas, sizeof(ParticulaConCodigo));
    orden.nuevo_orden = reservar_memoria(num_particulas, sizeof(int));
    orden.codigos = reservar_memoria(num_particulas, sizeof(uint64_t));
    orden.auxiliar = crear_particulas(num_particulas);
    return orden;
}

void liberar_orden_morton(OrdenMorton *orden)
{
    free(orden->pares);
    free(orden->auxiliar_pares);
    free(orden->nuevo_orden);
    free(orden->codigos);
    liberar_particulas(&orden->auxiliar);
}

/* El cubo más pequeño que contiene todas las partículas. Cada hilo busca el
 * mínimo y el máximo de su parte y OpenMP combina los resultados
 * (reduction), así ningún hilo escribe en una variable compartida. */
void calcular_cubo_raiz(const Particulas *particulas, Vector3 *centro, double *mitad_lado)
{
    double min_x = INFINITY, min_y = INFINITY, min_z = INFINITY;
    double max_x = -INFINITY, max_y = -INFINITY, max_z = -INFINITY;

    #pragma omp parallel for reduction(min : min_x, min_y, min_z) reduction(max : max_x, max_y, max_z)
    for (int i = 0; i < particulas->cantidad; i++) {
        Vector3 posicion = particulas->posiciones[i];
        min_x = fmin(min_x, posicion.x);
        min_y = fmin(min_y, posicion.y);
        min_z = fmin(min_z, posicion.z);
        max_x = fmax(max_x, posicion.x);
        max_y = fmax(max_y, posicion.y);
        max_z = fmax(max_z, posicion.z);
    }

    *centro = (Vector3){(min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2};
    *mitad_lado = fmax(max_x - min_x, fmax(max_y - min_y, max_z - min_z)) / 2;
}

/* Rebanada (de 0 a 2^21 - 1) en la que cae la coordenada dentro del cubo. */
static uint32_t rebanada(double coordenada, double esquina, double escala)
{
    double posicion = (coordenada - esquina) * escala;
    uint32_t ultima = (1u << BITS_POR_EJE) - 1;
    /* La partícula que está justo en el borde superior daría 2^21. */
    return posicion >= ultima ? ultima : (uint32_t)posicion;
}

/* Intercala los bits: ... z1 y1 x1 z0 y0 x0. Cada grupo de tres bits es un
 * octante, empezando por la raíz en los bits más altos. El orden x, y, z
 * coincide con el número de hijo: bit 0 = x, bit 1 = y, bit 2 = z. */
static uint64_t intercalar_bits(uint32_t x, uint32_t y, uint32_t z)
{
    uint64_t codigo = 0;
    for (int b = 0; b < BITS_POR_EJE; b++) {
        codigo |= (uint64_t)((x >> b) & 1) << (3 * b);
        codigo |= (uint64_t)((y >> b) & 1) << (3 * b + 1);
        codigo |= (uint64_t)((z >> b) & 1) << (3 * b + 2);
    }
    return codigo;
}

static void calcular_codigos(const Particulas *particulas, Vector3 centro, double mitad_lado,
                             ParticulaConCodigo *pares)
{
    Vector3 esquina = {centro.x - mitad_lado, centro.y - mitad_lado, centro.z - mitad_lado};
    /* Si todas las partículas están en el mismo punto, el cubo mide 0. */
    double escala = mitad_lado > 0.0 ? (1u << BITS_POR_EJE) / (2.0 * mitad_lado) : 0.0;

    #pragma omp parallel for
    for (int i = 0; i < particulas->cantidad; i++) {
        Vector3 posicion = particulas->posiciones[i];
        pares[i].codigo = intercalar_bits(rebanada(posicion.x, esquina.x, escala),
                                          rebanada(posicion.y, esquina.y, escala),
                                          rebanada(posicion.z, esquina.z, escala));
        pares[i].indice = i;
    }
}

/* Con el índice como desempate, el orden es siempre el mismo aunque dos
 * partículas tengan el mismo código: el árbol no depende del número de hilos. */
static int comparar_pares(const void *a, const void *b)
{
    const ParticulaConCodigo *primero = a;
    const ParticulaConCodigo *segundo = b;
    if (primero->codigo != segundo->codigo) {
        return primero->codigo < segundo->codigo ? -1 : 1;
    }
    return (primero->indice > segundo->indice) - (primero->indice < segundo->indice);
}

/* Une dos mitades ya ordenadas: datos[0, mitad) y datos[mitad, cantidad). */
static void mezclar(ParticulaConCodigo *datos, int mitad, int cantidad, ParticulaConCodigo *auxiliar)
{
    int izquierda = 0, derecha = mitad, destino = 0;
    while (izquierda < mitad && derecha < cantidad) {
        if (comparar_pares(&datos[izquierda], &datos[derecha]) <= 0) {
            auxiliar[destino++] = datos[izquierda++];
        } else {
            auxiliar[destino++] = datos[derecha++];
        }
    }
    while (izquierda < mitad) auxiliar[destino++] = datos[izquierda++];
    while (derecha < cantidad) auxiliar[destino++] = datos[derecha++];
    memcpy(datos, auxiliar, cantidad * sizeof(ParticulaConCodigo));
}

/* Merge sort con tareas: las dos mitades se ordenan en paralelo y luego se
 * mezclan. El taskwait es la barrera: no se puede mezclar hasta que las dos
 * mitades estén listas. */
static void ordenar_mezclando(ParticulaConCodigo *datos, ParticulaConCodigo *auxiliar, int cantidad)
{
    if (cantidad <= UMBRAL_ORDENAMIENTO_SECUENCIAL) {
        qsort(datos, cantidad, sizeof(ParticulaConCodigo), comparar_pares);
        return;
    }
    int mitad = cantidad / 2;
    #pragma omp task
    ordenar_mezclando(datos, auxiliar, mitad);
    #pragma omp task
    ordenar_mezclando(datos + mitad, auxiliar + mitad, cantidad - mitad);
    #pragma omp taskwait
    mezclar(datos, mitad, cantidad, auxiliar);
}

void ordenar_por_morton(Particulas *particulas, Vector3 centro, double mitad_lado, OrdenMorton *orden)
{
    int cantidad = particulas->cantidad;
    calcular_codigos(particulas, centro, mitad_lado, orden->pares);

    /* Un hilo arranca la recursión y los demás toman las tareas que va creando. */
    #pragma omp parallel
    #pragma omp single
    ordenar_mezclando(orden->pares, orden->auxiliar_pares, cantidad);

    #pragma omp parallel for
    for (int i = 0; i < cantidad; i++) {
        orden->nuevo_orden[i] = orden->pares[i].indice;
        orden->codigos[i] = orden->pares[i].codigo;
    }
    reordenar_particulas(particulas, &orden->auxiliar, orden->nuevo_orden);
}
