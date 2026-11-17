#ifndef CRITERIO_H
#define CRITERIO_H

#include <stdbool.h>

#include "arbol.h"

/* Distancia de un valor al intervalo [minimo, maximo] (0 si está adentro).
 * A lo sumo una de las dos diferencias es positiva. Se escribe con "a > b ?
 * a : b" y no con fmax ni con if: fmax tiene que tratar los NaN de una forma
 * especial y gcc no la convierte en una instrucción simple, y los if generan
 * saltos que el procesador predice mal. Con fmax el cálculo de fuerzas era
 * 1,8 veces más lento. */
static inline double distancia_a_intervalo(double valor, double minimo, double maximo)
{
    double por_debajo = minimo - valor;
    double por_encima = valor - maximo;
    double fuera = por_debajo > por_encima ? por_debajo : por_encima;
    return fuera > 0.0 ? fuera : 0.0;
}

/* CRITERIO DE APERTURA: aquí se decide si una celda se trata como una sola
 * partícula en su centro de masa (true) o se abre (false).
 *
 * Enunciado: tamaño / distancia < θ, comparado con cuadrados para no sacar
 * raíces. La distancia se mide desde el centro de masa de la celda hasta el
 * punto más cercano de una caja:
 *   - al calcular fuerzas, la caja es una partícula (mínimo = máximo);
 *   - al armar el árbol esencial, la caja es todo el dominio de otro proceso.
 * Una partícula nunca está más cerca de la celda que la caja que la
 * contiene, así que toda celda que se envía resumida cumple también el
 * criterio para cada partícula del otro proceso. Si se cambia el criterio
 * aquí, los dos usos cambian juntos y esa garantía se mantiene.
 *
 * Va en un .h como static inline porque se evalúa millones de veces por
 * paso: así el compilador lo copia en cada lugar donde se usa, sin el
 * costo de una llamada. */
static inline bool celda_suficientemente_lejana(const Nodo *nodo, Vector3 caja_minima, Vector3 caja_maxima,
                                                double theta)
{
    Vector3 centro = nodo->centro_de_masa;
    double dx = distancia_a_intervalo(centro.x, caja_minima.x, caja_maxima.x);
    double dy = distancia_a_intervalo(centro.y, caja_minima.y, caja_maxima.y);
    double dz = distancia_a_intervalo(centro.z, caja_minima.z, caja_maxima.z);
    double distancia2 = dx * dx + dy * dy + dz * dz;
    double tamano = 2.0 * nodo->mitad_lado;
    return tamano * tamano < theta * theta * distancia2;
}

#endif
