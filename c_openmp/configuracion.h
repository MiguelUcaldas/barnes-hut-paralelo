#ifndef CONFIGURACION_H
#define CONFIGURACION_H

typedef struct {
    const char *entrada;
    char distribucion[64];             /* para el CSV; se deduce del nombre del archivo */
    int pasos;
    double dt;
    double theta;
    double epsilon;
    int cada_energia;                  /* 0 = solo al inicio y al final */
    const char *guardar_aceleraciones; /* NULL = no guardar */
    const char *guardar_final;         /* NULL = no guardar */
    const char *csv;
} Configuracion;

Configuracion leer_configuracion(int argc, char **argv);

#endif
