#ifndef PARTICULAS_H
#define PARTICULAS_H

typedef struct {
    double x, y, z;
} Vector3;

typedef struct {
    int cantidad;
    Vector3 *posiciones;
    Vector3 *velocidades;
    Vector3 *aceleraciones;
    double *potenciales;   /* φ_i = -Σ_j m_j / √(r² + ε²); la energía potencial es ½ Σ m_i φ_i */
    double *masas;
} Particulas;

Particulas leer_particulas(const char *ruta);
void escribir_particulas(const char *ruta, const Particulas *particulas);
void escribir_aceleraciones(const char *ruta, const Particulas *particulas);
void liberar_particulas(Particulas *particulas);

#endif
