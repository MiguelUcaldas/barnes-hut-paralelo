"""Verifica la correctitud de cualquiera de las cuatro versiones.

fuerzas: compara las aceleraciones que calculó una versión sobre un archivo
         de partículas con las de la suma directa sobre ese mismo archivo.
         Error de cada partícula: |a_BH - a_directa| / |a_directa|.
energia: calcula E = K + U exacta (suma directa) en el estado inicial y en
         el estado final que dejó una versión, y reporta
         |E_final - E_inicial| / |E_inicial|.

Calcular la energía aquí, y no confiar en la que imprime cada versión,
permite detectar también un error en el cálculo de energía de la versión.

Uso:
    python comun/verificar.py fuerzas --particulas datos/uniforme_2000.txt --aceleraciones acel_openmp.txt --epsilon 0.01
    python comun/verificar.py energia --inicial datos/uniforme_2000.txt --final final_openmp.txt --epsilon 0.01
"""

import argparse
import sys

import numpy as np

from formato import leer_aceleraciones, leer_particulas
from referencia_directa import aceleraciones_directas, energia_total

# Con θ = 0,5 Barnes-Hut da un error promedio cercano a 0,5 %. Un promedio
# por encima del 1 % con θ ≤ 0,5 indica un error en el código, no en θ.
TOLERANCIA_FUERZAS = 1e-2
TOLERANCIA_ENERGIA = 1e-3


def errores_relativos(aceleraciones, referencia):
    diferencias = np.linalg.norm(aceleraciones - referencia, axis=1)
    magnitudes = np.linalg.norm(referencia, axis=1)
    # Si la aceleración exacta es cero (por ejemplo, todas las partículas en
    # un mismo punto), se usa el error absoluto para no dividir 0/0.
    return diferencias / np.where(magnitudes > 0, magnitudes, 1.0)


def verificar_fuerzas(ruta_particulas, ruta_aceleraciones, epsilon, tolerancia):
    posiciones, _, masas = leer_particulas(ruta_particulas)
    aceleraciones = leer_aceleraciones(ruta_aceleraciones)
    if len(aceleraciones) != len(masas):
        print(f"Error: {len(masas)} partículas pero {len(aceleraciones)} aceleraciones")
        return False

    referencia = aceleraciones_directas(posiciones, masas, epsilon)
    errores = errores_relativos(aceleraciones, referencia)
    peor = int(np.argmax(errores))

    print(f"Fuerzas: N = {len(masas)}, ε = {epsilon}")
    print(f"  error relativo promedio: {errores.mean():.3e}")
    print(f"  error relativo mediano:  {np.median(errores):.3e}")
    print(f"  error relativo máximo:   {errores[peor]:.3e}  (partícula {peor})")
    aprobado = errores.mean() < tolerancia
    print(f"  {'CORRECTO' if aprobado else 'INCORRECTO'}: promedio {'<' if aprobado else '>='} {tolerancia:g}")
    return aprobado


def verificar_energia(ruta_inicial, ruta_final, epsilon, tolerancia):
    posiciones_ini, velocidades_ini, masas_ini = leer_particulas(ruta_inicial)
    posiciones_fin, velocidades_fin, masas_fin = leer_particulas(ruta_final)
    if not np.array_equal(masas_ini, masas_fin):
        print("Error: las masas del estado final no coinciden con las del inicial")
        return False

    energia_inicial = energia_total(posiciones_ini, velocidades_ini, masas_ini, epsilon)
    energia_final = energia_total(posiciones_fin, velocidades_fin, masas_fin, epsilon)
    error = abs(energia_final - energia_inicial) / abs(energia_inicial)

    print(f"Energía: N = {len(masas_ini)}, ε = {epsilon}")
    print(f"  E inicial: {energia_inicial:.10f}")
    print(f"  E final:   {energia_final:.10f}")
    print(f"  error relativo: {error:.3e}")
    aprobado = error < tolerancia
    print(f"  {'CORRECTO' if aprobado else 'INCORRECTO'}: error {'<' if aprobado else '>='} {tolerancia:g}")
    return aprobado


def main():
    parser = argparse.ArgumentParser(description="Verifica fuerzas y energía contra la suma directa.")
    subcomandos = parser.add_subparsers(dest="comando", required=True)

    comando_fuerzas = subcomandos.add_parser("fuerzas")
    comando_fuerzas.add_argument("--particulas", required=True, help="archivo de entrada que usó la versión")
    comando_fuerzas.add_argument("--aceleraciones", required=True, help="aceleraciones que escribió la versión")
    comando_fuerzas.add_argument("--tolerancia", type=float, default=TOLERANCIA_FUERZAS)

    comando_energia = subcomandos.add_parser("energia")
    comando_energia.add_argument("--inicial", required=True, help="archivo de entrada de la simulación")
    comando_energia.add_argument("--final", required=True, help="estado final que escribió la versión")
    comando_energia.add_argument("--tolerancia", type=float, default=TOLERANCIA_ENERGIA)

    for comando in (comando_fuerzas, comando_energia):
        comando.add_argument("--epsilon", type=float, required=True)
    args = parser.parse_args()

    if args.comando == "fuerzas":
        aprobado = verificar_fuerzas(args.particulas, args.aceleraciones, args.epsilon, args.tolerancia)
    else:
        aprobado = verificar_energia(args.inicial, args.final, args.epsilon, args.tolerancia)
    # El código de salida permite usar la verificación desde los scripts de bash.
    sys.exit(0 if aprobado else 1)


if __name__ == "__main__":
    main()
