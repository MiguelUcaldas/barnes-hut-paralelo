#ifndef CONFIGURACION_H
#define CONFIGURACION_H

/* Cómo se reparten las partículas entre hilos al calcular las fuerzas. */
typedef enum {
    REPARTO_ESTATICO,   /* la misma cantidad de partículas por hilo */
    REPARTO_DINAMICO,   /* bloques pequeños que los hilos toman a medida que terminan */
    REPARTO_COSTO,      /* la misma cantidad de interacciones (medidas en el paso anterior) por hilo */
} Reparto;

typedef struct {
    const char *entrada;
    char distribucion[64];             /* para el CSV; se deduce del nombre del archivo */
    int pasos;
    double dt;
    double theta;
    double epsilon;
    int hilos;
    Reparto reparto;
    int cada_energia;                  /* 0 = solo al inicio y al final */
    const char *guardar_aceleraciones; /* NULL = no guardar */
    const char *guardar_final;         /* NULL = no guardar */
    const char *csv;
} Configuracion;

Configuracion leer_configuracion(int argc, char **argv);
const char *nombre_del_reparto(Reparto reparto);

#endif
