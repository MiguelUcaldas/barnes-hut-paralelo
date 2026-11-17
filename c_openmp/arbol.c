#include "arbol.h"

#include <omp.h>
#include <stdlib.h>

#include "memoria.h"
#include "morton.h"

/* Un subárbol con más partículas que esto se construye en una tarea aparte.
 * Más pequeño no vale la pena: crear la tarea costaría más que construirlo. */
#define UMBRAL_TAREA 1024

Arbol crear_arbol(int num_particulas)
{
    Arbol arbol;
    /* Con partículas repartidas hay unos 3N nodos, más lo que quede sin usar
     * al final de los bloques de cada hilo. Si no alcanza, construir_arbol
     * duplica la capacidad y vuelve a empezar. */
    arbol.capacidad = 4 * num_particulas + 8;
    arbol.nodos = reservar_memoria(arbol.capacidad, sizeof(Nodo));
    arbol.num_nodos = 0;
    arbol.nodos_en_uso = 0;
    arbol.num_bloques = omp_get_max_threads();
    arbol.bloques = reservar_memoria(arbol.num_bloques, sizeof(BloqueDeNodos));
    return arbol;
}

void liberar_arbol(Arbol *arbol)
{
    free(arbol->nodos);
    free(arbol->bloques);
}

static void iniciar_nodo(Nodo *nodo, Vector3 centro, double mitad_lado)
{
    *nodo = (Nodo){
        .centro = centro,
        .mitad_lado = mitad_lado,
        .masa = 0.0,
        .centro_de_masa = centro,
        .primer_hijo = SIN_HIJOS,
        .primera_particula = 0,
        .num_particulas = 0,
    };
}

/* Reserva 8 nodos seguidos para los hijos.
 *
 * Varios hilos construyen a la vez y no pueden recibir el mismo lugar. Si
 * cada reserva fuera un incremento atómico del contador global, todos los
 * hilos pelearían por esa variable (su línea de caché viajaría de núcleo en
 * núcleo en cada subdivisión) y la construcción no escalaría. Por eso cada
 * hilo toma del contador un bloque de NODOS_POR_BLOQUE nodos con una sola
 * operación atómica y reparte de ahí sin sincronizar: su bloque es solo suyo.
 *
 * Las tareas de OpenMP son "tied": entre leer omp_get_thread_num() y usar el
 * bloque no hay ningún punto donde el hilo pueda cambiar de tarea.
 *
 * Si ya no hay espacio se devuelve SIN_HIJOS; el contador queda por encima
 * de la capacidad y construir_arbol amplía el arreglo y repite. */
static int reservar_hijos(Arbol *arbol)
{
    BloqueDeNodos *bloque = &arbol->bloques[omp_get_thread_num()];
    if (bloque->siguiente == bloque->fin) {
        int inicio;
        #pragma omp atomic capture
        {
            inicio = arbol->num_nodos;
            arbol->num_nodos += NODOS_POR_BLOQUE;
        }
        bloque->siguiente = inicio;
        bloque->fin = inicio + NODOS_POR_BLOQUE;
    }
    int primer_hijo = bloque->siguiente;
    bloque->siguiente += 8;
    return primer_hijo + 8 <= arbol->capacidad ? primer_hijo : SIN_HIJOS;
}

static void calcular_masa_de_hoja(Nodo *nodo, const Particulas *particulas)
{
    double masa = 0.0;
    Vector3 momento = {0.0, 0.0, 0.0};   /* Σ m·x, para dividir al final */
    int ultima = nodo->primera_particula + nodo->num_particulas;
    for (int j = nodo->primera_particula; j < ultima; j++) {
        double m = particulas->masas[j];
        masa += m;
        momento.x += m * particulas->posiciones[j].x;
        momento.y += m * particulas->posiciones[j].y;
        momento.z += m * particulas->posiciones[j].z;
    }
    nodo->masa = masa;
    if (masa > 0.0) {
        nodo->centro_de_masa = (Vector3){momento.x / masa, momento.y / masa, momento.z / masa};
    }
}

static void sumar_masas_de_hijos(Nodo *nodo, const Nodo *nodos)
{
    double masa = 0.0;
    Vector3 momento = {0.0, 0.0, 0.0};
    for (int k = 0; k < 8; k++) {
        const Nodo *hijo = &nodos[nodo->primer_hijo + k];
        masa += hijo->masa;
        momento.x += hijo->masa * hijo->centro_de_masa.x;
        momento.y += hijo->masa * hijo->centro_de_masa.y;
        momento.z += hijo->masa * hijo->centro_de_masa.z;
    }
    nodo->masa = masa;
    if (masa > 0.0) {
        nodo->centro_de_masa = (Vector3){momento.x / masa, momento.y / masa, momento.z / masa};
    }
}

/* Los tres bits del código que dicen en qué hijo de un nodo de este nivel
 * cae la partícula. La raíz (nivel 0) usa los tres bits más altos. */
static int octante_en_nivel(uint64_t codigo, int nivel)
{
    return (int)((codigo >> (3 * (BITS_POR_EJE - 1 - nivel))) & 7);
}

/* Arma el subárbol del nodo con las partículas [primera, primera + cantidad).
 * Como están ordenadas por código, las del hijo 0 van primero, luego las del
 * hijo 1, etc.: basta recorrer el rango para cortarlo en 8 pedazos.
 * Cada subárbol escribe solo en sus propios nodos, así que dos tareas nunca
 * escriben en el mismo lugar; lo único compartido es el contador de nodos. */
static void construir_nodo(Arbol *arbol, const Particulas *particulas, const uint64_t *codigos,
                           int indice, int primera, int cantidad, int nivel)
{
    Nodo *nodo = &arbol->nodos[indice];
    nodo->primera_particula = primera;
    nodo->num_particulas = cantidad;

    if (cantidad <= 1 || nivel == BITS_POR_EJE) {
        calcular_masa_de_hoja(nodo, particulas);
        return;
    }
    int primer_hijo = reservar_hijos(arbol);
    if (primer_hijo == SIN_HIJOS) {
        return;   /* sin espacio: construir_arbol lo detecta y repite */
    }
    nodo->primer_hijo = primer_hijo;

    double cuarto = nodo->mitad_lado / 2;
    int siguiente = primera;
    int fin = primera + cantidad;
    for (int k = 0; k < 8; k++) {
        int desde = siguiente;
        while (siguiente < fin && octante_en_nivel(codigos[siguiente], nivel) == k) {
            siguiente++;
        }
        Vector3 centro_hijo = {
            nodo->centro.x + ((k & 1) ? cuarto : -cuarto),
            nodo->centro.y + ((k & 2) ? cuarto : -cuarto),
            nodo->centro.z + ((k & 4) ? cuarto : -cuarto),
        };
        iniciar_nodo(&arbol->nodos[primer_hijo + k], centro_hijo, cuarto);

        int cantidad_hijo = siguiente - desde;
        if (cantidad_hijo > UMBRAL_TAREA) {
            #pragma omp task
            construir_nodo(arbol, particulas, codigos, primer_hijo + k, desde, cantidad_hijo, nivel + 1);
        } else {
            construir_nodo(arbol, particulas, codigos, primer_hijo + k, desde, cantidad_hijo, nivel + 1);
        }
    }
    /* La masa del nodo necesita la de sus 8 hijos: hay que esperar a que
     * terminen las tareas creadas arriba. */
    #pragma omp taskwait
    sumar_masas_de_hijos(nodo, arbol->nodos);
}

void construir_arbol(Arbol *arbol, const Particulas *particulas, const uint64_t *codigos,
                     Vector3 centro, double mitad_lado)
{
    while (1) {
        arbol->num_nodos = 1;
        iniciar_nodo(&arbol->nodos[0], centro, mitad_lado);
        for (int hilo = 0; hilo < arbol->num_bloques; hilo++) {
            arbol->bloques[hilo].siguiente = arbol->bloques[hilo].fin = 0;
        }

        /* Un hilo empieza en la raíz; los demás toman las tareas de los
         * subárboles grandes a medida que aparecen. Un cúmulo denso genera
         * más tareas, y las toma el hilo que esté libre: así se balancea. */
        #pragma omp parallel
        #pragma omp single
        construir_nodo(arbol, particulas, codigos, 0, 0, particulas->cantidad, 0);

        if (arbol->num_nodos <= arbol->capacidad) {
            arbol->nodos_en_uso = arbol->num_nodos;
            for (int hilo = 0; hilo < arbol->num_bloques; hilo++) {
                arbol->nodos_en_uso -= arbol->bloques[hilo].fin - arbol->bloques[hilo].siguiente;
            }
            return;
        }
        /* Faltó espacio. No se puede usar realloc mientras otros hilos
         * escriben en el arreglo, así que se amplía aquí, ya sin hilos
         * trabajando, y se construye de nuevo. La capacidad queda para los
         * pasos siguientes, así que esto pasa muy pocas veces. */
        while (arbol->capacidad < arbol->num_nodos) {
            arbol->capacidad *= 2;
        }
        arbol->nodos = ampliar_memoria(arbol->nodos, arbol->capacidad, sizeof(Nodo));
    }
}
