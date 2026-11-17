#ifndef PARTICULAS_H
#define PARTICULAS_H

typedef struct {
    double x, y, z;
} Vector3;

/* En cada paso las partículas se reacomodan en orden de Morton, así que la
 * posición i del arreglo no es la fila i del archivo: esa está en
 * identificadores[i]. */
typedef struct {
    int cantidad;
    Vector3 *posiciones;
    Vector3 *velocidades;
    Vector3 *aceleraciones;
    double *potenciales;    /* φ_i = -Σ_j m_j / √(r² + ε²); la energía potencial es ½ Σ m_i φ_i */
    double *masas;
    int *interacciones;     /* cuántas atracciones sumó la partícula en el último cálculo: su costo */
    int *identificadores;   /* fila del archivo de entrada */
} Particulas;

Particulas crear_particulas(int cantidad);
Particulas leer_particulas(const char *ruta);
void reordenar_particulas(Particulas *particulas, Particulas *auxiliar, const int *nuevo_orden);
void escribir_particulas(const char *ruta, const Particulas *particulas);
void escribir_aceleraciones(const char *ruta, const Particulas *particulas);
void liberar_particulas(Particulas *particulas);

#endif
