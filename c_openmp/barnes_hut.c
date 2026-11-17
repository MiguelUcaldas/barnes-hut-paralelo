/* Simulación de N cuerpos con Barnes-Hut en C con OpenMP.
 *
 * Cada paso es un leapfrog patada-deriva-patada:
 *   1. v ← v + a·Δt/2
 *   2. x ← x + v·Δt
 *   3. ordenar por Morton, reconstruir el árbol y recalcular a
 *   4. v ← v + a·Δt/2
 *
 * Uso: bin/barnes_hut --entrada datos/uniforme_2000.txt --hilos 6 --reparto costo
 */

#include <omp.h>
#include <stdio.h>

#include "arbol.h"
#include "configuracion.h"
#include "fuerzas.h"
#include "medicion.h"
#include "morton.h"
#include "particulas.h"

/* Todo lo que se reutiliza de un paso al siguiente. */
typedef struct {
    Particulas particulas;
    Arbol arbol;
    OrdenMorton orden;
    RepartoDeTrabajo reparto;
} Simulacion;

static void media_patada(Particulas *particulas, double dt)
{
    #pragma omp parallel for
    for (int i = 0; i < particulas->cantidad; i++) {
        particulas->velocidades[i].x += particulas->aceleraciones[i].x * dt / 2;
        particulas->velocidades[i].y += particulas->aceleraciones[i].y * dt / 2;
        particulas->velocidades[i].z += particulas->aceleraciones[i].z * dt / 2;
    }
}

static void deriva(Particulas *particulas, double dt)
{
    #pragma omp parallel for
    for (int i = 0; i < particulas->cantidad; i++) {
        particulas->posiciones[i].x += particulas->velocidades[i].x * dt;
        particulas->posiciones[i].y += particulas->velocidades[i].y * dt;
        particulas->posiciones[i].z += particulas->velocidades[i].z * dt;
    }
}

/* Paso 3 del leapfrog: "foto congelada" de las posiciones actuales. Cada
 * fase es una región paralela distinta y termina con una barrera implícita:
 * el árbol no se arma hasta que el orden esté listo, y las fuerzas no se
 * calculan hasta que el árbol esté completo. */
static void recalcular_aceleraciones(Simulacion *sim, const Configuracion *config, TiemposPorFase *tiempos)
{
    Vector3 centro;
    double mitad_lado;

    double marca = segundos_actuales();
    calcular_cubo_raiz(&sim->particulas, &centro, &mitad_lado);
    tiempos->construccion += segundos_actuales() - marca;

    marca = segundos_actuales();
    ordenar_por_morton(&sim->particulas, centro, mitad_lado, &sim->orden);
    tiempos->ordenamiento += segundos_actuales() - marca;

    marca = segundos_actuales();
    construir_arbol(&sim->arbol, &sim->particulas, sim->orden.codigos, centro, mitad_lado);
    tiempos->construccion += segundos_actuales() - marca;

    marca = segundos_actuales();
    tiempos->desbalance = calcular_aceleraciones(&sim->arbol, &sim->particulas, config->theta, config->epsilon,
                                                 &sim->reparto);
    tiempos->fuerzas += segundos_actuales() - marca;
}

static TiemposPorFase avanzar_un_paso(Simulacion *sim, const Configuracion *config)
{
    TiemposPorFase tiempos = {0};
    double inicio_paso = segundos_actuales();

    double marca = segundos_actuales();
    media_patada(&sim->particulas, config->dt);
    deriva(&sim->particulas, config->dt);
    tiempos.integracion = segundos_actuales() - marca;

    recalcular_aceleraciones(sim, config, &tiempos);

    marca = segundos_actuales();
    media_patada(&sim->particulas, config->dt);
    tiempos.integracion += segundos_actuales() - marca;

    tiempos.total = segundos_actuales() - inicio_paso;
    return tiempos;
}

int main(int argc, char **argv)
{
    Configuracion config = leer_configuracion(argc, argv);
    /* Sin ajuste dinámico, cada región paralela tiene exactamente
     * config.hilos hilos; el reparto estático y por costo cuenta con eso. */
    omp_set_dynamic(0);
    omp_set_num_threads(config.hilos);

    Simulacion sim;
    sim.particulas = leer_particulas(config.entrada);
    sim.arbol = crear_arbol(sim.particulas.cantidad);
    sim.orden = crear_orden_morton(sim.particulas.cantidad);
    sim.reparto = crear_reparto(config.reparto, config.hilos);

    char version[64];
    snprintf(version, sizeof version, "c_openmp_%s", nombre_del_reparto(config.reparto));
    imprimir_entorno(version, &config, sim.particulas.cantidad);

    /* La primera media patada necesita las aceleraciones del estado inicial.
     * Después, las del final de un paso sirven para el comienzo del siguiente. */
    TiemposPorFase tiempos_iniciales = {0};
    recalcular_aceleraciones(&sim, &config, &tiempos_iniciales);
    printf("Cálculo inicial de fuerzas: %.3f s (%d nodos)\n\n",
           tiempos_iniciales.construccion + tiempos_iniciales.ordenamiento + tiempos_iniciales.fuerzas,
           sim.arbol.nodos_en_uso);
    if (config.guardar_aceleraciones != NULL) {
        escribir_aceleraciones(config.guardar_aceleraciones, &sim.particulas);
    }
    double energia_inicial = energia_total(&sim.particulas);

    TiemposPorFase acumulado = {0};
    imprimir_encabezado_pasos();
    for (int paso = 1; paso <= config.pasos; paso++) {
        TiemposPorFase tiempos = avanzar_un_paso(&sim, &config);
        sumar_tiempos(&acumulado, tiempos);
        imprimir_paso(paso, tiempos);
        if (config.cada_energia > 0 && paso % config.cada_energia == 0) {
            imprimir_energia(paso, energia_total(&sim.particulas), energia_inicial);
        }
    }

    double energia_final = energia_total(&sim.particulas);
    imprimir_resumen(acumulado, config.pasos, energia_inicial, energia_final);
    printf("  árbol del último paso: %d nodos\n", sim.arbol.nodos_en_uso);
    agregar_fila_csv(&config, version, sim.particulas.cantidad, acumulado, energia_inicial, energia_final);

    if (config.guardar_final != NULL) {
        escribir_particulas(config.guardar_final, &sim.particulas);
    }
    liberar_reparto(&sim.reparto);
    liberar_orden_morton(&sim.orden);
    liberar_arbol(&sim.arbol);
    liberar_particulas(&sim.particulas);
    return 0;
}
