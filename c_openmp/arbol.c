#include "arbol.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "memoria.h"

Arbol crear_arbol(int num_particulas)
{
    Arbol arbol;
    /* Con partículas repartidas hay unos 2N a 3N nodos. Si no alcanza,
     * nuevo_nodo duplica la capacidad. */
    arbol.capacidad = 4 * num_particulas + 8;
    arbol.nodos = reservar_memoria(arbol.capacidad, sizeof(Nodo));
    arbol.siguiente = reservar_memoria(num_particulas, sizeof(int));
    arbol.num_nodos = 0;
    arbol.profundidad_maxima = 0;
    return arbol;
}

void liberar_arbol(Arbol *arbol)
{
    free(arbol->nodos);
    free(arbol->siguiente);
}

/* El cubo más pequeño que contiene todas las partículas. Se recalcula en cada
 * paso porque las partículas se mueven. */
static void calcular_cubo_raiz(const Particulas *particulas, Vector3 *centro, double *mitad_lado)
{
    Vector3 minimo = particulas->posiciones[0];
    Vector3 maximo = minimo;
    for (int i = 1; i < particulas->cantidad; i++) {
        Vector3 posicion = particulas->posiciones[i];
        minimo.x = fmin(minimo.x, posicion.x);
        minimo.y = fmin(minimo.y, posicion.y);
        minimo.z = fmin(minimo.z, posicion.z);
        maximo.x = fmax(maximo.x, posicion.x);
        maximo.y = fmax(maximo.y, posicion.y);
        maximo.z = fmax(maximo.z, posicion.z);
    }
    *centro = (Vector3){(minimo.x + maximo.x) / 2, (minimo.y + maximo.y) / 2, (minimo.z + maximo.z) / 2};
    double lado = fmax(maximo.x - minimo.x, fmax(maximo.y - minimo.y, maximo.z - minimo.z));
    *mitad_lado = lado / 2;
}

static int nuevo_nodo(Arbol *arbol, Vector3 centro, double mitad_lado)
{
    if (arbol->num_nodos == arbol->capacidad) {
        arbol->capacidad *= 2;
        arbol->nodos = ampliar_memoria(arbol->nodos, arbol->capacidad, sizeof(Nodo));
    }
    int indice = arbol->num_nodos++;
    arbol->nodos[indice] = (Nodo){
        .centro = centro,
        .mitad_lado = mitad_lado,
        .masa = 0.0,
        .centro_de_masa = centro,
        .primer_hijo = SIN_HIJOS,
        .particula = SIN_PARTICULA,
    };
    return indice;
}

/* Número de octante de 0 a 7: un bit por eje (1 = lado positivo).
 * El hijo k tiene su centro desplazado hacia el lado que indican sus bits. */
static int octante(Vector3 posicion, Vector3 centro)
{
    return (posicion.x >= centro.x) | (posicion.y >= centro.y) << 1 | (posicion.z >= centro.z) << 2;
}

static void subdividir(Arbol *arbol, int indice)
{
    /* nuevo_nodo puede mover el arreglo con realloc: por eso se copian los
     * datos del padre antes y después se usa el índice, nunca un puntero
     * guardado de antes. */
    Vector3 centro = arbol->nodos[indice].centro;
    double cuarto = arbol->nodos[indice].mitad_lado / 2;
    int primer_hijo = arbol->num_nodos;
    for (int k = 0; k < 8; k++) {
        Vector3 centro_hijo = {
            centro.x + ((k & 1) ? cuarto : -cuarto),
            centro.y + ((k & 2) ? cuarto : -cuarto),
            centro.z + ((k & 4) ? cuarto : -cuarto),
        };
        nuevo_nodo(arbol, centro_hijo, cuarto);
    }
    arbol->nodos[indice].primer_hijo = primer_hijo;
}

/* Baja desde la raíz hasta una hoja. Si la hoja está vacía, la partícula se
 * queda ahí. Si está ocupada, la hoja se divide en 8, la partícula que estaba
 * baja a su octante y se sigue bajando con la nueva. */
static void insertar(Arbol *arbol, const Particulas *particulas, int i)
{
    Vector3 posicion = particulas->posiciones[i];
    int actual = 0;
    int profundidad = 0;
    while (true) {
        Nodo *nodo = &arbol->nodos[actual];
        if (nodo->primer_hijo != SIN_HIJOS) {
            actual = nodo->primer_hijo + octante(posicion, nodo->centro);
            profundidad++;
            continue;
        }
        if (nodo->particula == SIN_PARTICULA) {
            nodo->particula = i;
            break;
        }
        if (profundidad >= PROFUNDIDAD_MAXIMA) {
            /* Partículas en el mismo punto: dividir no las separaría nunca.
             * Se encadenan en esta hoja. */
            arbol->siguiente[i] = nodo->particula;
            nodo->particula = i;
            break;
        }
        int anterior = nodo->particula;
        nodo->particula = SIN_PARTICULA;
        subdividir(arbol, actual);
        Nodo *dividido = &arbol->nodos[actual];
        int hijo = dividido->primer_hijo + octante(particulas->posiciones[anterior], dividido->centro);
        arbol->nodos[hijo].particula = anterior;
    }
    if (profundidad > arbol->profundidad_maxima) {
        arbol->profundidad_maxima = profundidad;
    }
}

/* Los hijos siempre se crean después que su padre, así que tienen índices
 * mayores. Al recorrer el arreglo de atrás hacia adelante, cada nodo se
 * procesa después de todos sus hijos: es un recorrido de abajo hacia arriba
 * sin recursión. */
static void calcular_masas(Arbol *arbol, const Particulas *particulas)
{
    for (int indice = arbol->num_nodos - 1; indice >= 0; indice--) {
        Nodo *nodo = &arbol->nodos[indice];
        double masa = 0.0;
        Vector3 momento = {0.0, 0.0, 0.0};   /* Σ m·x, para dividir al final */

        if (nodo->primer_hijo == SIN_HIJOS) {
            for (int j = nodo->particula; j != SIN_PARTICULA; j = arbol->siguiente[j]) {
                double m = particulas->masas[j];
                masa += m;
                momento.x += m * particulas->posiciones[j].x;
                momento.y += m * particulas->posiciones[j].y;
                momento.z += m * particulas->posiciones[j].z;
            }
        } else {
            for (int k = 0; k < 8; k++) {
                const Nodo *hijo = &arbol->nodos[nodo->primer_hijo + k];
                masa += hijo->masa;
                momento.x += hijo->masa * hijo->centro_de_masa.x;
                momento.y += hijo->masa * hijo->centro_de_masa.y;
                momento.z += hijo->masa * hijo->centro_de_masa.z;
            }
        }

        nodo->masa = masa;
        if (masa > 0.0) {
            nodo->centro_de_masa = (Vector3){momento.x / masa, momento.y / masa, momento.z / masa};
        }
    }
}

void construir_arbol(Arbol *arbol, const Particulas *particulas)
{
    arbol->num_nodos = 0;
    arbol->profundidad_maxima = 0;
    for (int i = 0; i < particulas->cantidad; i++) {
        arbol->siguiente[i] = SIN_PARTICULA;
    }

    Vector3 centro;
    double mitad_lado;
    calcular_cubo_raiz(particulas, &centro, &mitad_lado);
    nuevo_nodo(arbol, centro, mitad_lado);

    for (int i = 0; i < particulas->cantidad; i++) {
        insertar(arbol, particulas, i);
    }
    calcular_masas(arbol, particulas);
}
