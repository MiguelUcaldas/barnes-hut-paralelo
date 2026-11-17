#include "fuerzas.h"

#include <math.h>
#include <stdbool.h>

#include "criterio.h"

#define NINGUNA (-1)

typedef struct {
    Vector3 aceleracion;
    double potencial;
    int interacciones;
} Acumulado;

/* Atracción de una masa puntual (partícula o celda) con suavizado:
 * a = m d / (|d|² + ε²)^(3/2) y φ = -m / √(|d|² + ε²). */
static void sumar_atraccion(Vector3 posicion, Vector3 fuente, double masa, double epsilon2, Acumulado *acumulado)
{
    double dx = fuente.x - posicion.x;
    double dy = fuente.y - posicion.y;
    double dz = fuente.z - posicion.z;
    double distancia2 = dx * dx + dy * dy + dz * dz + epsilon2;
    double inversa = 1.0 / sqrt(distancia2);
    double inversa3 = inversa * inversa * inversa;

    acumulado->aceleracion.x += masa * inversa3 * dx;
    acumulado->aceleracion.y += masa * inversa3 * dy;
    acumulado->aceleracion.z += masa * inversa3 * dz;
    acumulado->potencial -= masa * inversa;
    acumulado->interacciones++;
}

/* Sirve igual para el árbol local y para un árbol esencial recibido.
 * `excluida` es la fuente que es la propia partícula (solo en el local). */
static void recorrer_arbol(const Arbol *arbol, int indice_nodo, Vector3 posicion, int excluida,
                           double theta, double epsilon2, Acumulado *acumulado)
{
    const Nodo *nodo = &arbol->nodos[indice_nodo];
    if (nodo->masa == 0.0) {
        return;   /* celda vacía */
    }
    bool es_hoja = nodo->primer_hijo == SIN_HIJOS;
    /* Una hoja resumida viene de otro proceso, que ya comprobó que está lejos
     * de todo este dominio; no tiene hijos ni partículas que abrir. */
    bool resumida = es_hoja && nodo->num_particulas == 0;

    if (resumida || celda_suficientemente_lejana(nodo, posicion, posicion, theta)) {
        sumar_atraccion(posicion, nodo->centro_de_masa, nodo->masa, epsilon2, acumulado);
        return;
    }
    if (es_hoja) {
        for (int j = nodo->primera_particula; j < nodo->primera_particula + nodo->num_particulas; j++) {
            if (j != excluida) {
                sumar_atraccion(posicion, arbol->fuentes[j].posicion, arbol->fuentes[j].masa, epsilon2, acumulado);
            }
        }
        return;
    }
    for (int k = 0; k < 8; k++) {
        recorrer_arbol(arbol, nodo->primer_hijo + k, posicion, excluida, theta, epsilon2, acumulado);
    }
}

void calcular_aceleraciones(ParticulasLocales *particulas, const Arbol *local, const Arbol *remotos,
                            int num_procesos, int rango, double theta, double epsilon)
{
    double epsilon2 = epsilon * epsilon;
    for (int i = 0; i < particulas->cantidad; i++) {
        Particula *particula = &particulas->datos[i];
        Acumulado acumulado = {{0.0, 0.0, 0.0}, 0.0, 0};

        recorrer_arbol(local, 0, particula->posicion, i, theta, epsilon2, &acumulado);
        for (int r = 0; r < num_procesos; r++) {
            if (r != rango) {
                recorrer_arbol(&remotos[r], 0, particula->posicion, NINGUNA, theta, epsilon2, &acumulado);
            }
        }

        particula->aceleracion = acumulado.aceleracion;
        particula->potencial = acumulado.potencial;
        particula->costo = acumulado.interacciones;
    }
}

double energia_cinetica_local(const ParticulasLocales *particulas)
{
    double suma = 0.0;
    for (int i = 0; i < particulas->cantidad; i++) {
        Vector3 v = particulas->datos[i].velocidad;
        suma += particulas->datos[i].masa * (v.x * v.x + v.y * v.y + v.z * v.z);
    }
    return 0.5 * suma;
}

/* Cada par (i, j) aparece en φ_i y en φ_j, por eso el ½. */
double energia_potencial_local(const ParticulasLocales *particulas)
{
    double suma = 0.0;
    for (int i = 0; i < particulas->cantidad; i++) {
        suma += particulas->datos[i].masa * particulas->datos[i].potencial;
    }
    return 0.5 * suma;
}
