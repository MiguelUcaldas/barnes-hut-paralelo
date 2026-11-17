"""Suma directa O(N²): la referencia exacta contra la que se verifica Barnes-Hut.

Cada partícula recibe la atracción de todas las demás, sin aproximaciones.
Con G = 1 y suavizado ε:
    a_i = Σ_j m_j (x_j - x_i) / (|x_j - x_i|² + ε²)^(3/2)
    U   = -Σ_{i<j} m_i m_j / √(|x_j - x_i|² + ε²)
    K   = Σ_i m_i |v_i|² / 2

Uso:
    python comun/referencia_directa.py aceleraciones --particulas datos/uniforme_2000.txt --epsilon 0.01 --salida acel.txt
    python comun/referencia_directa.py energia --particulas datos/uniforme_2000.txt --epsilon 0.01
    python comun/referencia_directa.py simular --particulas datos/uniforme_2000.txt --epsilon 0.01 --dt 0.001 --pasos 100
"""

import argparse
import time

import numpy as np

from formato import escribir_aceleraciones, escribir_particulas, leer_particulas

# Los pares (i, j) se procesan por bloques de filas para que la memoria no
# crezca como N². Con 2 millones de pares, cada arreglo del bloque ocupa
# unos 50 MB, sin importar N.
PARES_POR_BLOQUE = 2_000_000


def bloques_de_filas(num_particulas):
    filas_por_bloque = max(1, PARES_POR_BLOQUE // num_particulas)
    for inicio in range(0, num_particulas, filas_por_bloque):
        yield inicio, min(inicio + filas_por_bloque, num_particulas)


def diferencias_del_bloque(posiciones, inicio, fin):
    """diferencias[b, j] = x_j - x_i con i = inicio + b: apunta de i hacia j."""
    return posiciones[None, :, :] - posiciones[inicio:fin, None, :]


def anular_autointeraccion(matriz, inicio, fin):
    # Una partícula no actúa sobre sí misma. Con ε = 0 ese término sería 0/0.
    filas = np.arange(fin - inicio)
    matriz[filas, inicio + filas] = 0.0


def aceleraciones_directas(posiciones, masas, epsilon):
    aceleraciones = np.zeros_like(posiciones)
    for inicio, fin in bloques_de_filas(len(masas)):
        diferencias = diferencias_del_bloque(posiciones, inicio, fin)
        distancias2 = np.einsum("bjk,bjk->bj", diferencias, diferencias) + epsilon**2
        with np.errstate(divide="ignore", invalid="ignore"):
            factores = masas[None, :] / distancias2**1.5
        anular_autointeraccion(factores, inicio, fin)
        aceleraciones[inicio:fin] = np.einsum("bj,bjk->bk", factores, diferencias)
    return aceleraciones


def energia_cinetica(velocidades, masas):
    return 0.5 * np.sum(masas * np.sum(velocidades**2, axis=1))


def energia_potencial(posiciones, masas, epsilon):
    suma = 0.0
    for inicio, fin in bloques_de_filas(len(masas)):
        diferencias = diferencias_del_bloque(posiciones, inicio, fin)
        distancias = np.sqrt(np.einsum("bjk,bjk->bj", diferencias, diferencias) + epsilon**2)
        with np.errstate(divide="ignore"):
            inversas = 1.0 / distancias
        anular_autointeraccion(inversas, inicio, fin)
        suma += masas[inicio:fin] @ inversas @ masas
    # Cada par (i, j) se contó dos veces: una desde i y otra desde j.
    return -0.5 * suma


def energia_total(posiciones, velocidades, masas, epsilon):
    return energia_cinetica(velocidades, masas) + energia_potencial(posiciones, masas, epsilon)


def simular(posiciones, velocidades, masas, epsilon, dt, pasos, cada=10):
    """Leapfrog patada-deriva-patada con fuerzas exactas.

    Sirve para saber cuánto error de energía pone el integrador solo, con
    un Δt dado; lo que una versión de Barnes-Hut tenga por encima de esto
    viene de la aproximación del árbol.
    """
    energia_inicial = energia_total(posiciones, velocidades, masas, epsilon)
    print(f"{'paso':>6} {'energía':>16} {'error relativo':>15}")
    print(f"{0:>6} {energia_inicial:>16.10f} {0.0:>15.3e}")

    # Las aceleraciones del estado inicial se calculan una vez antes del
    # ciclo: la primera media patada las necesita.
    aceleraciones = aceleraciones_directas(posiciones, masas, epsilon)
    for paso in range(1, pasos + 1):
        velocidades += aceleraciones * (dt / 2)
        posiciones += velocidades * dt
        aceleraciones = aceleraciones_directas(posiciones, masas, epsilon)
        velocidades += aceleraciones * (dt / 2)

        if paso % cada == 0 or paso == pasos:
            energia = energia_total(posiciones, velocidades, masas, epsilon)
            error = abs(energia - energia_inicial) / abs(energia_inicial)
            print(f"{paso:>6} {energia:>16.10f} {error:>15.3e}")
    return posiciones, velocidades


def main():
    parser = argparse.ArgumentParser(description="Suma directa O(N²) con G = 1 y suavizado ε.")
    subcomandos = parser.add_subparsers(dest="comando", required=True)

    comando_aceleraciones = subcomandos.add_parser("aceleraciones", help="escribe las aceleraciones exactas (ax ay az)")
    comando_energia = subcomandos.add_parser("energia", help="imprime K, U y E = K + U")
    comando_simular = subcomandos.add_parser("simular", help="leapfrog con fuerzas exactas")
    for comando in (comando_aceleraciones, comando_energia, comando_simular):
        comando.add_argument("--particulas", required=True)
        comando.add_argument("--epsilon", type=float, required=True)
    comando_aceleraciones.add_argument("--salida", required=True)
    comando_simular.add_argument("--dt", type=float, required=True)
    comando_simular.add_argument("--pasos", type=int, required=True)
    comando_simular.add_argument("--cada", type=int, default=10, help="imprimir la energía cada tantos pasos")
    comando_simular.add_argument("--salida", help="archivo para el estado final (opcional)")
    args = parser.parse_args()

    posiciones, velocidades, masas = leer_particulas(args.particulas)
    print(f"Archivo: {args.particulas}, N = {len(masas)}, ε = {args.epsilon}")

    if args.comando == "aceleraciones":
        inicio = time.perf_counter()
        aceleraciones = aceleraciones_directas(posiciones, masas, args.epsilon)
        print(f"Tiempo: {time.perf_counter() - inicio:.2f} s")
        escribir_aceleraciones(args.salida, aceleraciones)
        print(f"Aceleraciones en {args.salida}")

    elif args.comando == "energia":
        cinetica = energia_cinetica(velocidades, masas)
        potencial = energia_potencial(posiciones, masas, args.epsilon)
        print(f"K = {cinetica:.10f}")
        print(f"U = {potencial:.10f}")
        print(f"E = {cinetica + potencial:.10f}")
        print(f"2K/|U| = {2 * cinetica / abs(potencial):.4f}  (1 = equilibrio)")

    elif args.comando == "simular":
        print(f"Δt = {args.dt}, pasos = {args.pasos}")
        posiciones, velocidades = simular(posiciones, velocidades, masas, args.epsilon,
                                          args.dt, args.pasos, args.cada)
        if args.salida:
            escribir_particulas(args.salida, posiciones, velocidades, masas)
            print(f"Estado final en {args.salida}")


if __name__ == "__main__":
    main()
