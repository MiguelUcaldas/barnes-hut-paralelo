#ifndef MEDICION_H
#define MEDICION_H

#include "configuracion.h"

/* Segundos por fase. Cada proceso mide los suyos; para imprimir se toma el
 * máximo entre procesos (el más lento marca el ritmo). En el resumen son la
 * suma de todos los pasos. */
typedef struct {
    double ordenamiento;    /* códigos Morton, cortes y envío de partículas a su dueño */
    double construccion;    /* árbol local */
    double esencial;        /* extraer y enviar los árboles esenciales */
    double fuerzas;         /* recorrido, hasta que el último proceso termina */
    double integracion;
    double total;
    double desbalance;      /* proceso más lento / promedio, en el cálculo de fuerzas */
} TiemposPorFase;

double error_relativo(double valor, double referencia);
void sumar_tiempos(TiemposPorFase *acumulado, TiemposPorFase paso);

/* Operación colectiva. Devuelve, en el proceso 0, el máximo de cada fase
 * entre procesos y el desbalance de `calculo_de_fuerzas` (el tiempo de
 * cálculo propio, sin la espera en la barrera). */
TiemposPorFase reunir_tiempos(TiemposPorFase propios, double calculo_de_fuerzas);

void imprimir_entorno(const char *version, const Configuracion *config, int num_particulas, int num_procesos);
void imprimir_encabezado_pasos(void);
void imprimir_paso(int paso, TiemposPorFase tiempos);
void imprimir_energia(int paso, double energia, double energia_inicial);
void imprimir_resumen(TiemposPorFase acumulado, int pasos, double energia_inicial, double energia_final);
void agregar_fila_csv(const Configuracion *config, const char *version, int num_particulas, int num_procesos,
                      TiemposPorFase acumulado, double energia_inicial, double energia_final);

#endif
