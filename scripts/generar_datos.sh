#!/usr/bin/env bash
# Genera todos los archivos de partículas del proyecto con la semilla fija.
#   datos/<distribucion>_2000.txt          pequeños, para verificar contra la suma directa
#   datos/grandes/<distribucion>_<N>.txt   para medir (no se suben a git)
# Uso, desde la raíz del repositorio:
#   bash scripts/generar_datos.sh [ruta_de_python]
set -euo pipefail

PYTHON="${1:-.venv/bin/python}"
DISTRIBUCIONES="uniforme dos_cumulos concentrada"

for distribucion in $DISTRIBUCIONES; do
    "$PYTHON" comun/generar_particulas.py --distribucion "$distribucion" --n 2000
    for n in 20000 50000 100000 200000; do
        "$PYTHON" comun/generar_particulas.py --distribucion "$distribucion" --n "$n" \
            --salida "datos/grandes/${distribucion}_${n}.txt"
    done
done
