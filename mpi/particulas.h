#ifndef PARTICULAS_H
#define PARTICULAS_H

#include <stdint.h>

typedef struct {
    double x, y, z;
} Vector3;

/* Todo el estado de una partícula en un solo bloque de memoria, para poder
 * enviarla a otro proceso con un solo mensaje cuando cambia de dueño. */
typedef struct {
    Vector3 posicion;
    Vector3 velocidad;
    Vector3 aceleracion;
    double masa;
    double potencial;      /* φ_i; la energía potencial es ½ Σ m_i φ_i */
    double costo;          /* interacciones que sumó en el último cálculo de fuerzas */
    uint64_t codigo;       /* código Morton en el paso actual */
    int identificador;     /* fila del archivo de entrada */
} Particula;

/* Lo único que el árbol necesita de una partícula para atraer a otras.
 * Es también lo que se envía de las hojas del árbol esencial. */
typedef struct {
    Vector3 posicion;
    double masa;
} Fuente;

/* Las partículas de las que este proceso es dueño. Su cantidad cambia en
 * cada paso, cuando se reparten de nuevo. */
typedef struct {
    Particula *datos;
    int cantidad;
} ParticulasLocales;

#endif
