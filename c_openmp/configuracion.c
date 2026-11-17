#include "configuracion.h"

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void mostrar_ayuda(const char *programa)
{
    printf("Uso: %s --entrada archivo [opciones]\n\n"
           "  --entrada archivo           partículas: x y z vx vy vz masa (obligatorio)\n"
           "  --pasos N                   pasos de tiempo (100)\n"
           "  --dt valor                  paso de tiempo Δt (0.001)\n"
           "  --theta valor               criterio de apertura θ (0.5)\n"
           "  --epsilon valor             suavizado ε (0.01)\n"
           "  --cada-energia K            imprimir la energía cada K pasos (0 = solo inicio y final)\n"
           "  --guardar-aceleraciones f   aceleraciones del estado inicial, para verificar.py\n"
           "  --guardar-final f           estado final, para verificar.py\n"
           "  --distribucion nombre       nombre para el CSV (por defecto, sale del archivo)\n"
           "  --csv archivo               dónde agregar la medición (resultados/mediciones.csv)\n",
           programa);
}

/* "datos/grandes/dos_cumulos_20000.txt" -> "dos_cumulos": se quita la
 * carpeta, la extensión y el _N final. */
static void distribucion_desde_ruta(const char *ruta, char *destino, size_t tamano)
{
    const char *barra = strrchr(ruta, '/');
    snprintf(destino, tamano, "%s", barra ? barra + 1 : ruta);

    char *punto = strrchr(destino, '.');
    if (punto != NULL) {
        *punto = '\0';
    }
    char *guion = strrchr(destino, '_');
    if (guion != NULL && guion[1] != '\0' && strspn(guion + 1, "0123456789") == strlen(guion + 1)) {
        *guion = '\0';
    }
}

static void error_de_opcion(const char *mensaje)
{
    fprintf(stderr, "Error: %s (use --ayuda)\n", mensaje);
    exit(1);
}

Configuracion leer_configuracion(int argc, char **argv)
{
    Configuracion config = {
        .entrada = NULL,
        .distribucion = "",
        .pasos = 100,
        .dt = 0.001,
        .theta = 0.5,
        .epsilon = 0.01,
        .cada_energia = 0,
        .guardar_aceleraciones = NULL,
        .guardar_final = NULL,
        .csv = "resultados/mediciones.csv",
    };

    static const struct option opciones[] = {
        {"entrada", required_argument, NULL, 'e'},
        {"pasos", required_argument, NULL, 'p'},
        {"dt", required_argument, NULL, 'd'},
        {"theta", required_argument, NULL, 't'},
        {"epsilon", required_argument, NULL, 's'},
        {"cada-energia", required_argument, NULL, 'k'},
        {"guardar-aceleraciones", required_argument, NULL, 'a'},
        {"guardar-final", required_argument, NULL, 'f'},
        {"distribucion", required_argument, NULL, 'n'},
        {"csv", required_argument, NULL, 'c'},
        {"ayuda", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0},
    };

    int opcion;
    while ((opcion = getopt_long(argc, argv, "", opciones, NULL)) != -1) {
        switch (opcion) {
        case 'e': config.entrada = optarg; break;
        case 'p': config.pasos = atoi(optarg); break;
        case 'd': config.dt = atof(optarg); break;
        case 't': config.theta = atof(optarg); break;
        case 's': config.epsilon = atof(optarg); break;
        case 'k': config.cada_energia = atoi(optarg); break;
        case 'a': config.guardar_aceleraciones = optarg; break;
        case 'f': config.guardar_final = optarg; break;
        case 'n': snprintf(config.distribucion, sizeof config.distribucion, "%s", optarg); break;
        case 'c': config.csv = optarg; break;
        case 'h': mostrar_ayuda(argv[0]); exit(0);
        default: error_de_opcion("opción no reconocida");
        }
    }

    if (config.entrada == NULL) error_de_opcion("falta --entrada");
    if (config.pasos < 0) error_de_opcion("--pasos debe ser 0 o más");
    if (config.dt <= 0.0) error_de_opcion("--dt debe ser positivo");
    if (config.theta < 0.0) error_de_opcion("--theta no puede ser negativo");
    if (config.epsilon < 0.0) error_de_opcion("--epsilon no puede ser negativo");

    if (config.distribucion[0] == '\0') {
        distribucion_desde_ruta(config.entrada, config.distribucion, sizeof config.distribucion);
    }
    return config;
}
