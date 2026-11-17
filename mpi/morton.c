#include "morton.h"

/* Rebanada (de 0 a 2^21 - 1) en la que cae la coordenada dentro del cubo. */
static uint32_t rebanada(double coordenada, double esquina, double escala)
{
    double posicion = (coordenada - esquina) * escala;
    uint32_t ultima = (1u << BITS_POR_EJE) - 1;
    if (posicion <= 0.0) {
        return 0;
    }
    /* La partícula que está justo en el borde superior daría 2^21. */
    return posicion >= ultima ? ultima : (uint32_t)posicion;
}

/* Intercala los bits: ... z1 y1 x1 z0 y0 x0. Cada grupo de tres bits es un
 * octante, empezando por la raíz en los bits más altos; bit 0 = x, bit 1 = y,
 * bit 2 = z, igual que la numeración de los hijos. */
static uint64_t intercalar_bits(uint32_t x, uint32_t y, uint32_t z)
{
    uint64_t codigo = 0;
    for (int b = 0; b < BITS_POR_EJE; b++) {
        codigo |= (uint64_t)((x >> b) & 1) << (3 * b);
        codigo |= (uint64_t)((y >> b) & 1) << (3 * b + 1);
        codigo |= (uint64_t)((z >> b) & 1) << (3 * b + 2);
    }
    return codigo;
}

uint64_t codigo_morton(Vector3 posicion, Vector3 esquina, double escala)
{
    return intercalar_bits(rebanada(posicion.x, esquina.x, escala),
                           rebanada(posicion.y, esquina.y, escala),
                           rebanada(posicion.z, esquina.z, escala));
}

int octante_en_nivel(uint64_t codigo, int nivel)
{
    return (int)((codigo >> (3 * (BITS_POR_EJE - 1 - nivel))) & 7);
}
