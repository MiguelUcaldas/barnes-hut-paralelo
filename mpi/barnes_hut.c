/* Simulación de N cuerpos con Barnes-Hut en MPI.
 *
 * Cada proceso es dueño de un tramo contiguo de la curva de Morton (una
 * región compacta del espacio) y de las partículas que caen en él. En cada
 * paso de leapfrog patada-deriva-patada:
 *   1. v ← v + a·Δt/2 y 2. x ← x + v·Δt   (cada proceso, sus partículas)
 *   3. se reparten de nuevo las partículas por la curva de Morton, cada
 *      proceso construye el árbol de las suyas, recibe de los demás solo la
 *      parte de sus árboles que necesita (árbol esencial) y calcula a
 *   4. v ← v + a·Δt/2
 *
 * Uso: mpirun -np 6 mpi/bin/barnes_hut --entrada datos/uniforme_2000.txt --particion costo
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "arbol.h"
#include "configuracion.h"
#include "entrada_salida.h"
#include "esencial.h"
#include "fuerzas.h"
#include "medicion.h"
#include "memoria.h"
#include "particion.h"
#include "particulas.h"

/* Todo lo que se reutiliza de un paso al siguiente en este proceso. */
typedef struct {
    int rango;
    int num_procesos;
    int total;                 /* N global */
    ParticulasLocales particulas;
    Fuente *fuentes;           /* posición y masa de las partículas locales, para el árbol */
    uint64_t *codigos;
    int capacidad_fuentes;
    Arbol arbol;
    Intercambio intercambio;
} Simulacion;

static void media_patada(ParticulasLocales *particulas, double dt)
{
    for (int i = 0; i < particulas->cantidad; i++) {
        Particula *p = &particulas->datos[i];
        p->velocidad.x += p->aceleracion.x * dt / 2;
        p->velocidad.y += p->aceleracion.y * dt / 2;
        p->velocidad.z += p->aceleracion.z * dt / 2;
    }
}

static void deriva(ParticulasLocales *particulas, double dt)
{
    for (int i = 0; i < particulas->cantidad; i++) {
        Particula *p = &particulas->datos[i];
        p->posicion.x += p->velocidad.x * dt;
        p->posicion.y += p->velocidad.y * dt;
        p->posicion.z += p->velocidad.z * dt;
    }
}

/* El árbol solo necesita posición y masa; se copian a un arreglo compacto
 * para que el recorrido lea menos memoria. La fuente i es la partícula i. */
static void preparar_fuentes(Simulacion *sim)
{
    int cantidad = sim->particulas.cantidad;
    if (cantidad > sim->capacidad_fuentes) {
        sim->capacidad_fuentes = 2 * cantidad;
        sim->fuentes = ampliar_memoria(sim->fuentes, sim->capacidad_fuentes, sizeof(Fuente));
        sim->codigos = ampliar_memoria(sim->codigos, sim->capacidad_fuentes, sizeof(uint64_t));
    }
    for (int i = 0; i < cantidad; i++) {
        sim->fuentes[i] = (Fuente){sim->particulas.datos[i].posicion, sim->particulas.datos[i].masa};
        sim->codigos[i] = sim->particulas.datos[i].codigo;
    }
}

/* Paso 3 del leapfrog. Devuelve el tiempo de cálculo de fuerzas de este
 * proceso, sin la espera en la barrera, para medir el desbalance. */
static double recalcular_aceleraciones(Simulacion *sim, const Configuracion *config, TiemposPorFase *tiempos)
{
    Vector3 centro;
    double mitad_lado;

    double marca = MPI_Wtime();
    calcular_cubo_raiz_global(&sim->particulas, &centro, &mitad_lado);
    repartir_particulas(&sim->particulas, centro, mitad_lado, config->particion);
    tiempos->ordenamiento += MPI_Wtime() - marca;

    marca = MPI_Wtime();
    preparar_fuentes(sim);
    construir_arbol_local(&sim->arbol, sim->fuentes, sim->codigos, sim->particulas.cantidad, centro, mitad_lado);
    tiempos->construccion += MPI_Wtime() - marca;

    marca = MPI_Wtime();
    intercambiar_arboles_esenciales(&sim->intercambio, &sim->arbol, &sim->particulas, config->theta);
    tiempos->esencial += MPI_Wtime() - marca;

    marca = MPI_Wtime();
    calcular_aceleraciones(&sim->particulas, &sim->arbol, sim->intercambio.remotos, sim->num_procesos, sim->rango,
                           config->theta, config->epsilon);
    double calculo = MPI_Wtime() - marca;
    /* Barrera de la "foto congelada". En MPI no es necesaria para que el
     * resultado sea correcto: cada proceso solo mueve sus propias partículas,
     * y lo que recibió de los demás es una copia que nadie más modifica. Se
     * pone para que el tiempo de fuerzas sea el del proceso más lento y la
     * espera no se cuele en la fase siguiente. */
    MPI_Barrier(MPI_COMM_WORLD);
    tiempos->fuerzas += MPI_Wtime() - marca;
    return calculo;
}

static TiemposPorFase avanzar_un_paso(Simulacion *sim, const Configuracion *config, double *calculo)
{
    TiemposPorFase tiempos = {0};
    double inicio_paso = MPI_Wtime();

    double marca = MPI_Wtime();
    media_patada(&sim->particulas, config->dt);
    deriva(&sim->particulas, config->dt);
    tiempos.integracion = MPI_Wtime() - marca;

    *calculo = recalcular_aceleraciones(sim, config, &tiempos);

    marca = MPI_Wtime();
    media_patada(&sim->particulas, config->dt);
    tiempos.integracion += MPI_Wtime() - marca;

    tiempos.total = MPI_Wtime() - inicio_paso;
    return tiempos;
}

static double energia_total(const Simulacion *sim)
{
    double propia = energia_cinetica_local(&sim->particulas) + energia_potencial_local(&sim->particulas);
    double total;
    MPI_Allreduce(&propia, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return total;
}

/* Cuántas partículas tiene cada proceso y cuánto árbol recibió: muestra el
 * balance del particionado y que nadie recibe el árbol completo. */
static void imprimir_reparto(const Simulacion *sim)
{
    int cantidad = sim->particulas.cantidad;
    int minimo, maximo;
    MPI_Reduce(&cantidad, &minimo, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&cantidad, &maximo, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    double propios[3] = {sim->arbol.num_nodos, sim->intercambio.nodos_recibidos_total,
                         sim->intercambio.fuentes_recibidas_total};
    double sumas[3], maximos[3];
    MPI_Reduce(propios, sumas, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(propios, maximos, 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (sim->rango == 0) {
        int p = sim->num_procesos;
        printf("\nReparto en el último paso\n");
        printf("  partículas por proceso:        mínimo %d, máximo %d (promedio %d)\n", minimo, maximo, sim->total / p);
        printf("  nodos del árbol local:         promedio %.0f, total entre procesos %.0f\n", sumas[0] / p, sumas[0]);
        printf("  árbol esencial recibido:       promedio %.0f nodos y %.0f partículas (máximo %.0f nodos)\n",
               sumas[1] / p, sumas[2] / p, maximos[1]);
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    Simulacion sim;
    MPI_Comm_rank(MPI_COMM_WORLD, &sim.rango);
    MPI_Comm_size(MPI_COMM_WORLD, &sim.num_procesos);

    Configuracion config = leer_configuracion(argc, argv, sim.rango);
    sim.particulas = leer_y_repartir(config.entrada, &sim.total);
    sim.capacidad_fuentes = 0;
    sim.fuentes = NULL;
    sim.codigos = NULL;
    sim.arbol = crear_arbol(4 * (sim.total / sim.num_procesos) + 8);
    sim.intercambio = crear_intercambio();

    char version[64];
    snprintf(version, sizeof version, "mpi_%s", nombre_de_particion(config.particion));
    if (sim.rango == 0) {
        imprimir_entorno(version, &config, sim.total, sim.num_procesos);
    }

    /* La primera media patada necesita las aceleraciones del estado inicial. */
    TiemposPorFase iniciales = {0};
    double inicio = MPI_Wtime();
    recalcular_aceleraciones(&sim, &config, &iniciales);
    if (sim.rango == 0) {
        printf("Cálculo inicial de fuerzas: %.3f s\n\n", MPI_Wtime() - inicio);
    }
    if (config.guardar_aceleraciones != NULL) {
        escribir_aceleraciones_globales(config.guardar_aceleraciones, &sim.particulas, sim.total);
    }
    double energia_inicial = energia_total(&sim);

    TiemposPorFase acumulado = {0};
    if (sim.rango == 0) {
        imprimir_encabezado_pasos();
    }
    for (int paso = 1; paso <= config.pasos; paso++) {
        double calculo;
        TiemposPorFase propios = avanzar_un_paso(&sim, &config, &calculo);
        TiemposPorFase tiempos = reunir_tiempos(propios, calculo);
        if (sim.rango == 0) {
            sumar_tiempos(&acumulado, tiempos);
            imprimir_paso(paso, tiempos);
        }
        if (config.cada_energia > 0 && paso % config.cada_energia == 0) {
            double energia = energia_total(&sim);
            if (sim.rango == 0) {
                imprimir_energia(paso, energia, energia_inicial);
            }
        }
    }

    double energia_final = energia_total(&sim);
    if (sim.rango == 0) {
        imprimir_resumen(acumulado, config.pasos, energia_inicial, energia_final);
    }
    imprimir_reparto(&sim);
    if (sim.rango == 0) {
        agregar_fila_csv(&config, version, sim.total, sim.num_procesos, acumulado, energia_inicial, energia_final);
    }
    if (config.guardar_final != NULL) {
        escribir_estado_global(config.guardar_final, &sim.particulas, sim.total);
    }

    liberar_intercambio(&sim.intercambio);
    liberar_arbol(&sim.arbol);
    free(sim.fuentes);
    free(sim.codigos);
    free(sim.particulas.datos);
    MPI_Finalize();
    return 0;
}
