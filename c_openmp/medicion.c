#include "medicion.h"

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

/* CLOCK_MONOTONIC no salta si el reloj del sistema se ajusta durante la
 * corrida, a diferencia de la hora del día. */
double segundos_actuales(void)
{
    struct timespec ahora;
    clock_gettime(CLOCK_MONOTONIC, &ahora);
    return ahora.tv_sec + ahora.tv_nsec * 1e-9;
}

double error_relativo(double valor, double referencia)
{
    return fabs(valor - referencia) / fabs(referencia);
}

void sumar_tiempos(TiemposPorFase *acumulado, TiemposPorFase paso)
{
    acumulado->construccion += paso.construccion;
    acumulado->ordenamiento += paso.ordenamiento;
    acumulado->fuerzas += paso.fuerzas;
    acumulado->integracion += paso.integracion;
    acumulado->total += paso.total;
}

void imprimir_entorno(const char *version, const Configuracion *config, int num_particulas, int trabajadores)
{
    printf("Barnes-Hut, versión %s\n", version);
    printf("  entrada:             %s (distribución %s)\n", config->entrada, config->distribucion);
    printf("  N:                   %d\n", num_particulas);
    printf("  pasos:               %d\n", config->pasos);
    printf("  θ:                   %g\n", config->theta);
    printf("  ε:                   %g\n", config->epsilon);
    printf("  Δt:                  %g\n", config->dt);
    printf("  trabajadores:        %d\n", trabajadores);
    printf("  núcleos disponibles: %ld\n", sysconf(_SC_NPROCESSORS_ONLN));
    printf("\n");
}

/* printf cuenta bytes y no letras: en UTF-8 "ó" ocupa 2 bytes. Para que las
 * columnas con tildes queden alineadas, se suma al ancho un byte por cada
 * byte de continuación (los que empiezan con los bits 10). */
static int bytes_de_mas(const char *texto)
{
    int extra = 0;
    for (const unsigned char *c = (const unsigned char *)texto; *c != '\0'; c++) {
        if ((*c & 0xC0) == 0x80) {
            extra++;
        }
    }
    return extra;
}

void imprimir_encabezado_pasos(void)
{
    const char *columnas[] = {"paso", "construcción", "ordenamiento", "fuerzas", "integración", "total"};
    const int anchos[] = {6, 13, 13, 10, 12, 10};
    printf("Tiempos por paso (ms)\n");
    for (int k = 0; k < 6; k++) {
        printf("%s%*s", k > 0 ? " " : "", anchos[k] + bytes_de_mas(columnas[k]), columnas[k]);
    }
    printf("\n");
}

void imprimir_paso(int paso, TiemposPorFase tiempos)
{
    printf("%6d %13.2f %13.2f %10.2f %12.2f %10.2f\n", paso,
           1e3 * tiempos.construccion, 1e3 * tiempos.ordenamiento, 1e3 * tiempos.fuerzas,
           1e3 * tiempos.integracion, 1e3 * tiempos.total);
}

void imprimir_energia(int paso, double energia, double energia_inicial)
{
    printf("       energía en el paso %d: %.10f (error relativo %.3e)\n",
           paso, energia, error_relativo(energia, energia_inicial));
}

static void imprimir_fase(const char *nombre, double segundos, int pasos, double total)
{
    double por_paso = pasos > 0 ? 1e3 * segundos / pasos : 0.0;
    double porcentaje = total > 0.0 ? 100.0 * segundos / total : 0.0;
    printf("  %-*s %10.3f %14.2f %11.1f %%\n", 24 + bytes_de_mas(nombre), nombre, segundos, por_paso, porcentaje);
}

void imprimir_resumen(TiemposPorFase acumulado, int pasos, double energia_inicial, double energia_final)
{
    printf("\nTiempos acumulados de %d pasos\n", pasos);
    printf("  %-24s %10s %14s %13s\n", "fase", "total (s)", "por paso (ms)", "porcentaje");
    imprimir_fase("construcción del árbol", acumulado.construccion, pasos, acumulado.total);
    imprimir_fase("ordenamiento", acumulado.ordenamiento, pasos, acumulado.total);
    imprimir_fase("fuerzas", acumulado.fuerzas, pasos, acumulado.total);
    imprimir_fase("integración", acumulado.integracion, pasos, acumulado.total);
    imprimir_fase("total", acumulado.total, pasos, acumulado.total);

    printf("\nEnergía (calculada con el árbol)\n");
    printf("  inicial:        %.10f\n", energia_inicial);
    printf("  final:          %.10f\n", energia_final);
    printf("  error relativo: %.3e\n", error_relativo(energia_final, energia_inicial));
}

void agregar_fila_csv(const Configuracion *config, const char *version, int num_particulas, int trabajadores,
                      TiemposPorFase acumulado, double energia_inicial, double energia_final)
{
    FILE *archivo = fopen(config->csv, "a");
    if (archivo == NULL) {
        fprintf(stderr, "Aviso: no se pudo abrir %s; la medición no se guardó\n", config->csv);
        return;
    }
    /* El encabezado solo va si el archivo está vacío, así varias corridas
     * de distintas versiones quedan en la misma tabla. */
    fseek(archivo, 0, SEEK_END);
    if (ftell(archivo) == 0) {
        fprintf(archivo, "version,distribucion,N,pasos,theta,epsilon,dt,trabajadores,"
                         "t_construccion,t_ordenamiento,t_fuerzas,t_integracion,t_total,"
                         "energia_inicial,energia_final,error_energia\n");
    }
    fprintf(archivo, "%s,%s,%d,%d,%g,%g,%g,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.10g,%.10g,%.6e\n",
            version, config->distribucion, num_particulas, config->pasos,
            config->theta, config->epsilon, config->dt, trabajadores,
            acumulado.construccion, acumulado.ordenamiento, acumulado.fuerzas,
            acumulado.integracion, acumulado.total,
            energia_inicial, energia_final, error_relativo(energia_final, energia_inicial));
    fclose(archivo);
    printf("\nMedición agregada a %s\n", config->csv);
}
