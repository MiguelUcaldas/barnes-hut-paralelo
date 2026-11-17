#include "arbol.h"

#include <stdlib.h>

#include "memoria.h"
#include "morton.h"

Arbol crear_arbol(int capacidad_inicial)
{
    Arbol arbol;
    arbol.capacidad = capacidad_inicial < 8 ? 8 : capacidad_inicial;
    arbol.nodos = reservar_memoria(arbol.capacidad, sizeof(Nodo));
    arbol.num_nodos = 0;
    arbol.fuentes = NULL;
    return arbol;
}

void liberar_arbol(Arbol *arbol)
{
    free(arbol->nodos);
}

static void iniciar_nodo(Nodo *nodo, Vector3 centro, double mitad_lado)
{
    *nodo = (Nodo){
        .centro = centro,
        .mitad_lado = mitad_lado,
        .masa = 0.0,
        .centro_de_masa = centro,
        .primer_hijo = SIN_HIJOS,
        .primera_particula = 0,
        .num_particulas = 0,
    };
}

/* Reserva 8 nodos seguidos. Puede mover el arreglo con realloc, así que
 * después de llamarla solo se usan índices, nunca punteros guardados antes. */
static int reservar_hijos(Arbol *arbol)
{
    if (arbol->num_nodos + 8 > arbol->capacidad) {
        arbol->capacidad *= 2;
        arbol->nodos = ampliar_memoria(arbol->nodos, arbol->capacidad, sizeof(Nodo));
    }
    int primer_hijo = arbol->num_nodos;
    arbol->num_nodos += 8;
    return primer_hijo;
}

static void calcular_masa_de_hoja(Nodo *nodo, const Fuente *fuentes)
{
    double masa = 0.0;
    Vector3 momento = {0.0, 0.0, 0.0};   /* Σ m·x, para dividir al final */
    for (int j = nodo->primera_particula; j < nodo->primera_particula + nodo->num_particulas; j++) {
        masa += fuentes[j].masa;
        momento.x += fuentes[j].masa * fuentes[j].posicion.x;
        momento.y += fuentes[j].masa * fuentes[j].posicion.y;
        momento.z += fuentes[j].masa * fuentes[j].posicion.z;
    }
    nodo->masa = masa;
    if (masa > 0.0) {
        nodo->centro_de_masa = (Vector3){momento.x / masa, momento.y / masa, momento.z / masa};
    }
}

static void sumar_masas_de_hijos(Nodo *nodo, const Nodo *nodos)
{
    double masa = 0.0;
    Vector3 momento = {0.0, 0.0, 0.0};
    for (int k = 0; k < 8; k++) {
        const Nodo *hijo = &nodos[nodo->primer_hijo + k];
        masa += hijo->masa;
        momento.x += hijo->masa * hijo->centro_de_masa.x;
        momento.y += hijo->masa * hijo->centro_de_masa.y;
        momento.z += hijo->masa * hijo->centro_de_masa.z;
    }
    nodo->masa = masa;
    if (masa > 0.0) {
        nodo->centro_de_masa = (Vector3){momento.x / masa, momento.y / masa, momento.z / masa};
    }
}

/* Igual que en la versión OpenMP, sin tareas: las fuentes están ordenadas
 * por código, así que las de cada hijo son un tramo contiguo del rango. */
static void construir_nodo(Arbol *arbol, const uint64_t *codigos, int indice, int primera, int cantidad, int nivel)
{
    arbol->nodos[indice].primera_particula = primera;
    arbol->nodos[indice].num_particulas = cantidad;
    if (cantidad <= 1 || nivel == BITS_POR_EJE) {
        calcular_masa_de_hoja(&arbol->nodos[indice], arbol->fuentes);
        return;
    }

    int primer_hijo = reservar_hijos(arbol);
    arbol->nodos[indice].primer_hijo = primer_hijo;
    Vector3 centro = arbol->nodos[indice].centro;
    double cuarto = arbol->nodos[indice].mitad_lado / 2;

    int siguiente = primera;
    int fin = primera + cantidad;
    for (int k = 0; k < 8; k++) {
        int desde = siguiente;
        while (siguiente < fin && octante_en_nivel(codigos[siguiente], nivel) == k) {
            siguiente++;
        }
        Vector3 centro_hijo = {
            centro.x + ((k & 1) ? cuarto : -cuarto),
            centro.y + ((k & 2) ? cuarto : -cuarto),
            centro.z + ((k & 4) ? cuarto : -cuarto),
        };
        iniciar_nodo(&arbol->nodos[primer_hijo + k], centro_hijo, cuarto);
        construir_nodo(arbol, codigos, primer_hijo + k, desde, siguiente - desde, nivel + 1);
    }
    sumar_masas_de_hijos(&arbol->nodos[indice], arbol->nodos);
}

void construir_arbol_local(Arbol *arbol, const Fuente *fuentes, const uint64_t *codigos, int cantidad,
                           Vector3 centro, double mitad_lado)
{
    arbol->fuentes = fuentes;
    arbol->num_nodos = 1;
    iniciar_nodo(&arbol->nodos[0], centro, mitad_lado);
    construir_nodo(arbol, codigos, 0, 0, cantidad, 0);
}
