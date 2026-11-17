#ifndef MEDICION_H
#define MEDICION_H

#include "configuracion.h"

/* Segundos por fase. En un paso son los de ese paso; en el resumen, la suma
 * de todos los pasos. */
typedef struct {
    double construccion;
    double ordenamiento;
    double fuerzas;
    double integracion;
    double total;
} TiemposPorFase;

double segundos_actuales(void);
double error_relativo(double valor, double referencia);
void sumar_tiempos(TiemposPorFase *acumulado, TiemposPorFase paso);

void imprimir_entorno(const char *version, const Configuracion *config, int num_particulas, int trabajadores);
void imprimir_encabezado_pasos(void);
void imprimir_paso(int paso, TiemposPorFase tiempos);
void imprimir_energia(int paso, double energia, double energia_inicial);
void imprimir_resumen(TiemposPorFase acumulado, int pasos, double energia_inicial, double energia_final);
void agregar_fila_csv(const Configuracion *config, const char *version, int num_particulas, int trabajadores,
                      TiemposPorFase acumulado, double energia_inicial, double energia_final);

#endif
