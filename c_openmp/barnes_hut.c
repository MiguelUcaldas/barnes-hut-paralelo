/* Simulación de N cuerpos con Barnes-Hut, versión secuencial en C.
 *
 * Cada paso es un leapfrog patada-deriva-patada:
 *   1. v ← v + a·Δt/2
 *   2. x ← x + v·Δt
 *   3. reconstruir el árbol y recalcular a
 *   4. v ← v + a·Δt/2
 *
 * Uso: bin/barnes_hut --entrada datos/uniforme_2000.txt --pasos 100 --theta 0.5
 */

#include <stdio.h>

#include "arbol.h"
#include "configuracion.h"
#include "fuerzas.h"
#include "medicion.h"
#include "particulas.h"

#define VERSION "c_secuencial"
#define TRABAJADORES 1

static void media_patada(Particulas *particulas, double dt)
{
    for (int i = 0; i < particulas->cantidad; i++) {
        particulas->velocidades[i].x += particulas->aceleraciones[i].x * dt / 2;
        particulas->velocidades[i].y += particulas->aceleraciones[i].y * dt / 2;
        particulas->velocidades[i].z += particulas->aceleraciones[i].z * dt / 2;
    }
}

static void deriva(Particulas *particulas, double dt)
{
    for (int i = 0; i < particulas->cantidad; i++) {
        particulas->posiciones[i].x += particulas->velocidades[i].x * dt;
        particulas->posiciones[i].y += particulas->velocidades[i].y * dt;
        particulas->posiciones[i].z += particulas->velocidades[i].z * dt;
    }
}

static TiemposPorFase avanzar_un_paso(Arbol *arbol, Particulas *particulas, const Configuracion *config)
{
    TiemposPorFase tiempos = {0};
    double inicio_paso = segundos_actuales();

    double marca = segundos_actuales();
    media_patada(particulas, config->dt);
    deriva(particulas, config->dt);
    tiempos.integracion = segundos_actuales() - marca;

    /* "Foto congelada": el árbol se arma con las posiciones de este instante
     * y todas las aceleraciones se calculan con él antes de mover a nadie. */
    marca = segundos_actuales();
    construir_arbol(arbol, particulas);
    tiempos.construccion = segundos_actuales() - marca;

    marca = segundos_actuales();
    calcular_aceleraciones(arbol, particulas, config->theta, config->epsilon);
    tiempos.fuerzas = segundos_actuales() - marca;

    marca = segundos_actuales();
    media_patada(particulas, config->dt);
    tiempos.integracion += segundos_actuales() - marca;

    tiempos.total = segundos_actuales() - inicio_paso;
    return tiempos;
}

int main(int argc, char **argv)
{
    Configuracion config = leer_configuracion(argc, argv);
    Particulas particulas = leer_particulas(config.entrada);
    Arbol arbol = crear_arbol(particulas.cantidad);
    imprimir_entorno(VERSION, &config, particulas.cantidad, TRABAJADORES);

    /* La primera media patada necesita las aceleraciones del estado inicial.
     * Después, las del final de un paso sirven para el comienzo del siguiente. */
    double inicio = segundos_actuales();
    construir_arbol(&arbol, &particulas);
    calcular_aceleraciones(&arbol, &particulas, config.theta, config.epsilon);
    printf("Cálculo inicial de fuerzas: %.3f s (%d nodos, profundidad %d)\n\n",
           segundos_actuales() - inicio, arbol.num_nodos, arbol.profundidad_maxima);
    if (config.guardar_aceleraciones != NULL) {
        escribir_aceleraciones(config.guardar_aceleraciones, &particulas);
    }
    double energia_inicial = energia_total(&particulas);

    TiemposPorFase acumulado = {0};
    imprimir_encabezado_pasos();
    for (int paso = 1; paso <= config.pasos; paso++) {
        TiemposPorFase tiempos = avanzar_un_paso(&arbol, &particulas, &config);
        sumar_tiempos(&acumulado, tiempos);
        imprimir_paso(paso, tiempos);
        if (config.cada_energia > 0 && paso % config.cada_energia == 0) {
            imprimir_energia(paso, energia_total(&particulas), energia_inicial);
        }
    }

    double energia_final = energia_total(&particulas);
    imprimir_resumen(acumulado, config.pasos, energia_inicial, energia_final);
    printf("  árbol del último paso: %d nodos, profundidad %d\n", arbol.num_nodos, arbol.profundidad_maxima);
    agregar_fila_csv(&config, VERSION, particulas.cantidad, TRABAJADORES, acumulado, energia_inicial, energia_final);

    if (config.guardar_final != NULL) {
        escribir_particulas(config.guardar_final, &particulas);
    }
    liberar_arbol(&arbol);
    liberar_particulas(&particulas);
    return 0;
}
