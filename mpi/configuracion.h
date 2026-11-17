#ifndef CONFIGURACION_H
#define CONFIGURACION_H

/* Cómo se corta la curva de Morton entre procesos. */
typedef enum {
    PARTICION_VOLUMEN,    /* tramos de igual volumen (reparto ingenuo) */
    PARTICION_CANTIDAD,   /* tramos con el mismo número de partículas */
    PARTICION_COSTO,      /* tramos con el mismo número de interacciones (paso anterior) */
} Particion;

typedef struct {
    const char *entrada;
    char distribucion[64];             /* para el CSV; se deduce del nombre del archivo */
    int pasos;
    double dt;
    double theta;
    double epsilon;
    Particion particion;
    int cada_energia;                  /* 0 = solo al inicio y al final */
    const char *guardar_aceleraciones; /* NULL = no guardar */
    const char *guardar_final;         /* NULL = no guardar */
    const char *csv;
} Configuracion;

/* Todos los procesos leen las mismas opciones; solo el 0 imprime errores. */
Configuracion leer_configuracion(int argc, char **argv, int rango);
const char *nombre_de_particion(Particion particion);

#endif
