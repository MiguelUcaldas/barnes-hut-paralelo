#include "fuerzas.h"

#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdlib.h>

#include "memoria.h"

/* En el reparto dinámico, cada hilo toma bloques de este tamaño. Más pequeño
 * balancea mejor pero los hilos compiten más seguido por el siguiente bloque. */
#define BLOQUE_DINAMICO 64

typedef struct {
    Vector3 aceleracion;
    double potencial;
    int interacciones;
} Acumulado;

/* CRITERIO DE APERTURA: aquí se decide si una celda se abre o se trata como
 * una sola partícula en su centro de masa.
 * Enunciado: tamaño / distancia < θ. Se compara con cuadrados para no sacar
 * una raíz por cada celda visitada. Con θ = 0 nunca se cumple y el cálculo
 * baja hasta las hojas: es la suma directa. */
static bool celda_suficientemente_lejana(const Nodo *nodo, double distancia2, double theta)
{
    double tamano = 2.0 * nodo->mitad_lado;
    return tamano * tamano < theta * theta * distancia2;
}

/* Atracción de una masa puntual (partícula o celda) sobre la posición dada,
 * con suavizado: a = m d / (|d|² + ε²)^(3/2) y φ = -m / √(|d|² + ε²). */
static void sumar_atraccion(Vector3 posicion, Vector3 fuente, double masa, double epsilon2, Acumulado *acumulado)
{
    double dx = fuente.x - posicion.x;
    double dy = fuente.y - posicion.y;
    double dz = fuente.z - posicion.z;
    double distancia2 = dx * dx + dy * dy + dz * dz + epsilon2;
    double inversa = 1.0 / sqrt(distancia2);
    double inversa3 = inversa * inversa * inversa;

    acumulado->aceleracion.x += masa * inversa3 * dx;
    acumulado->aceleracion.y += masa * inversa3 * dy;
    acumulado->aceleracion.z += masa * inversa3 * dz;
    acumulado->potencial -= masa * inversa;
    acumulado->interacciones++;
}

static void recorrer_arbol(const Arbol *arbol, int indice_nodo, const Particulas *particulas, int i,
                           double theta, double epsilon2, Acumulado *acumulado)
{
    const Nodo *nodo = &arbol->nodos[indice_nodo];
    Vector3 posicion = particulas->posiciones[i];

    if (nodo->primer_hijo == SIN_HIJOS) {
        /* Hoja: interacción exacta con su partícula (o sus partículas, si
         * tienen el mismo código Morton), saltándose a sí misma. */
        int ultima = nodo->primera_particula + nodo->num_particulas;
        for (int j = nodo->primera_particula; j < ultima; j++) {
            if (j != i) {
                sumar_atraccion(posicion, particulas->posiciones[j], particulas->masas[j], epsilon2, acumulado);
            }
        }
        return;
    }

    double dx = nodo->centro_de_masa.x - posicion.x;
    double dy = nodo->centro_de_masa.y - posicion.y;
    double dz = nodo->centro_de_masa.z - posicion.z;
    double distancia2 = dx * dx + dy * dy + dz * dz;

    if (celda_suficientemente_lejana(nodo, distancia2, theta)) {
        sumar_atraccion(posicion, nodo->centro_de_masa, nodo->masa, epsilon2, acumulado);
        return;
    }
    for (int k = 0; k < 8; k++) {
        recorrer_arbol(arbol, nodo->primer_hijo + k, particulas, i, theta, epsilon2, acumulado);
    }
}

/* Solo lee el árbol y las posiciones, y solo escribe en la posición i: dos
 * hilos con distinto i nunca escriben en el mismo lugar. */
static void fuerza_sobre_particula(const Arbol *arbol, Particulas *particulas, int i, double theta, double epsilon2)
{
    Acumulado acumulado = {{0.0, 0.0, 0.0}, 0.0, 0};
    recorrer_arbol(arbol, 0, particulas, i, theta, epsilon2, &acumulado);
    particulas->aceleraciones[i] = acumulado.aceleracion;
    particulas->potenciales[i] = acumulado.potencial;
    particulas->interacciones[i] = acumulado.interacciones;
}

RepartoDeTrabajo crear_reparto(Reparto tipo, int num_hilos)
{
    RepartoDeTrabajo reparto;
    reparto.reparto = tipo;
    reparto.num_hilos = num_hilos;
    reparto.limites = reservar_memoria(num_hilos + 1, sizeof(int));
    reparto.tiempos = reservar_memoria(num_hilos, sizeof(double));
    return reparto;
}

void liberar_reparto(RepartoDeTrabajo *reparto)
{
    free(reparto->limites);
    free(reparto->tiempos);
}

/* Reparto por costo: rangos contiguos con la misma suma de interacciones.
 * El costo de cada partícula es el que tuvo en el paso anterior; como se
 * mueven poco entre pasos, es una buena predicción del costo actual. */
static void limites_por_costo(RepartoDeTrabajo *reparto, const Particulas *particulas)
{
    int n = particulas->cantidad;
    double costo_total = 0.0;
    for (int i = 0; i < n; i++) {
        costo_total += particulas->interacciones[i];
    }

    int hilo = 1;
    double acumulado = 0.0;
    reparto->limites[0] = 0;
    for (int i = 0; i < n && hilo < reparto->num_hilos; i++) {
        acumulado += particulas->interacciones[i];
        while (hilo < reparto->num_hilos && acumulado >= costo_total * hilo / reparto->num_hilos) {
            reparto->limites[hilo++] = i + 1;
        }
    }
    while (hilo <= reparto->num_hilos) {
        reparto->limites[hilo++] = n;
    }
}

/* REPARTO DE TRABAJO ENTRE HILOS: aquí se decide qué partículas calcula cada
 * hilo en los repartos estático y por costo (el dinámico lo hace OpenMP en
 * calcular_aceleraciones). Las partículas están en orden de Morton, así que
 * un rango contiguo es una región compacta del espacio: las partículas de un
 * hilo recorren casi los mismos nodos y los encuentran en la caché. */
static void calcular_limites(RepartoDeTrabajo *reparto, const Particulas *particulas)
{
    if (reparto->reparto == REPARTO_COSTO) {
        limites_por_costo(reparto, particulas);
        return;
    }
    for (int hilo = 0; hilo <= reparto->num_hilos; hilo++) {
        reparto->limites[hilo] = (int)((long)particulas->cantidad * hilo / reparto->num_hilos);
    }
}

static double desbalance(const RepartoDeTrabajo *reparto)
{
    double maximo = 0.0, suma = 0.0;
    for (int hilo = 0; hilo < reparto->num_hilos; hilo++) {
        maximo = fmax(maximo, reparto->tiempos[hilo]);
        suma += reparto->tiempos[hilo];
    }
    return suma > 0.0 ? maximo / (suma / reparto->num_hilos) : 1.0;
}

double calcular_aceleraciones(const Arbol *arbol, Particulas *particulas, double theta, double epsilon,
                              RepartoDeTrabajo *reparto)
{
    double epsilon2 = epsilon * epsilon;
    calcular_limites(reparto, particulas);

    /* El equipo tiene exactamente num_hilos hilos (omp_set_dynamic(0) en
     * main), así que cada hilo tiene su rango en limites. */
    #pragma omp parallel num_threads(reparto->num_hilos)
    {
        int hilo = omp_get_thread_num();
        double inicio = omp_get_wtime();

        if (reparto->reparto == REPARTO_DINAMICO) {
            #pragma omp for schedule(dynamic, BLOQUE_DINAMICO) nowait
            for (int i = 0; i < particulas->cantidad; i++) {
                fuerza_sobre_particula(arbol, particulas, i, theta, epsilon2);
            }
        } else {
            for (int i = reparto->limites[hilo]; i < reparto->limites[hilo + 1]; i++) {
                fuerza_sobre_particula(arbol, particulas, i, theta, epsilon2);
            }
        }
        reparto->tiempos[hilo] = omp_get_wtime() - inicio;
    }
    /* Al cerrar la región paralela hay una barrera implícita: ningún hilo
     * sigue (y ninguna partícula se mueve) hasta que todas las fuerzas estén
     * calculadas. Es la "foto congelada". */
    return desbalance(reparto);
}

double energia_cinetica(const Particulas *particulas)
{
    double suma = 0.0;
    #pragma omp parallel for reduction(+ : suma)
    for (int i = 0; i < particulas->cantidad; i++) {
        Vector3 v = particulas->velocidades[i];
        suma += particulas->masas[i] * (v.x * v.x + v.y * v.y + v.z * v.z);
    }
    return 0.5 * suma;
}

/* Usa los potenciales que dejó el último calcular_aceleraciones, así que no
 * cuesta otro recorrido del árbol. Cada par (i, j) aparece en φ_i y en φ_j,
 * por eso el ½. */
double energia_potencial(const Particulas *particulas)
{
    double suma = 0.0;
    #pragma omp parallel for reduction(+ : suma)
    for (int i = 0; i < particulas->cantidad; i++) {
        suma += particulas->masas[i] * particulas->potenciales[i];
    }
    return 0.5 * suma;
}

double energia_total(const Particulas *particulas)
{
    return energia_cinetica(particulas) + energia_potencial(particulas);
}
