#include "fuerzas.h"

#include <math.h>
#include <stdbool.h>

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
static void sumar_atraccion(Vector3 posicion, Vector3 fuente, double masa, double epsilon2,
                            Vector3 *aceleracion, double *potencial)
{
    double dx = fuente.x - posicion.x;
    double dy = fuente.y - posicion.y;
    double dz = fuente.z - posicion.z;
    double distancia2 = dx * dx + dy * dy + dz * dz + epsilon2;
    double inversa = 1.0 / sqrt(distancia2);
    double inversa3 = inversa * inversa * inversa;

    aceleracion->x += masa * inversa3 * dx;
    aceleracion->y += masa * inversa3 * dy;
    aceleracion->z += masa * inversa3 * dz;
    *potencial -= masa * inversa;
}

static void recorrer_arbol(const Arbol *arbol, int indice_nodo, const Particulas *particulas, int i,
                           double theta, double epsilon2, Vector3 *aceleracion, double *potencial)
{
    const Nodo *nodo = &arbol->nodos[indice_nodo];
    Vector3 posicion = particulas->posiciones[i];

    if (nodo->primer_hijo == SIN_HIJOS) {
        /* Hoja: interacción exacta con su partícula (o sus partículas, si
         * hay varias en el mismo punto), saltándose a sí misma. */
        for (int j = nodo->particula; j != SIN_PARTICULA; j = arbol->siguiente[j]) {
            if (j != i) {
                sumar_atraccion(posicion, particulas->posiciones[j], particulas->masas[j],
                                epsilon2, aceleracion, potencial);
            }
        }
        return;
    }

    double dx = nodo->centro_de_masa.x - posicion.x;
    double dy = nodo->centro_de_masa.y - posicion.y;
    double dz = nodo->centro_de_masa.z - posicion.z;
    double distancia2 = dx * dx + dy * dy + dz * dz;

    if (celda_suficientemente_lejana(nodo, distancia2, theta)) {
        sumar_atraccion(posicion, nodo->centro_de_masa, nodo->masa, epsilon2, aceleracion, potencial);
        return;
    }
    for (int k = 0; k < 8; k++) {
        recorrer_arbol(arbol, nodo->primer_hijo + k, particulas, i, theta, epsilon2, aceleracion, potencial);
    }
}

/* Cada partícula solo lee el árbol y las posiciones, y solo escribe su propia
 * aceleración: las iteraciones de este ciclo son independientes entre sí. */
void calcular_aceleraciones(const Arbol *arbol, Particulas *particulas, double theta, double epsilon)
{
    double epsilon2 = epsilon * epsilon;
    for (int i = 0; i < particulas->cantidad; i++) {
        Vector3 aceleracion = {0.0, 0.0, 0.0};
        double potencial = 0.0;
        recorrer_arbol(arbol, 0, particulas, i, theta, epsilon2, &aceleracion, &potencial);
        particulas->aceleraciones[i] = aceleracion;
        particulas->potenciales[i] = potencial;
    }
}

double energia_cinetica(const Particulas *particulas)
{
    double suma = 0.0;
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
    for (int i = 0; i < particulas->cantidad; i++) {
        suma += particulas->masas[i] * particulas->potenciales[i];
    }
    return 0.5 * suma;
}

double energia_total(const Particulas *particulas)
{
    return energia_cinetica(particulas) + energia_potencial(particulas);
}
