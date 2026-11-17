# Mapa del código

Dónde está cada concepto, para encontrarlo rápido durante la demostración. Se actualiza al terminar cada fase.

## Repositorio

| Concepto | Archivo | Función |
|---|---|---|
| Enunciado oficial | `docs/enunciado-proyecto-final.pdf` | — |
| Configuración de git (saltos de línea LF) | `.gitattributes` | — |
| Dependencias de Python | `requirements.txt` | — |
| Versiones y equipo | `ENTORNO.md` | — |

## Herramientas comunes (`comun/`)

| Concepto | Archivo | Función |
|---|---|---|
| Leer / escribir partículas (7 columnas) | `comun/formato.py` | `leer_particulas`, `escribir_particulas` |
| Leer / escribir aceleraciones (3 columnas) | `comun/formato.py` | `leer_aceleraciones`, `escribir_aceleraciones` |
| Distribución uniforme del enunciado | `comun/generar_particulas.py` | `uniforme` |
| Dos cúmulos en colisión | `comun/generar_particulas.py` | `dos_cumulos` (constantes `RADIO_CUMULO`, `SEPARACION_CUMULOS`, `VELOCIDAD_CHOQUE`) |
| Partículas concentradas | `comun/generar_particulas.py` | `concentrada` (constantes `FRACCION_CONCENTRADA`, `RADIO_CONCENTRADO`, `CENTRO_CONCENTRADO`) |
| Cúmulo en equilibrio (Plummer) | `comun/generar_particulas.py` | `esfera_de_plummer`, `fracciones_de_escape` |
| Barajar filas y masas 1/N | `comun/generar_particulas.py` | `generar` |
| Aceleraciones exactas O(N²) | `comun/referencia_directa.py` | `aceleraciones_directas` |
| Energía cinética / potencial / total | `comun/referencia_directa.py` | `energia_cinetica`, `energia_potencial`, `energia_total` |
| Leapfrog con fuerzas exactas | `comun/referencia_directa.py` | `simular` |
| Tamaño de los bloques de la suma directa | `comun/referencia_directa.py` | `PARES_POR_BLOQUE` |
| Error relativo de fuerzas | `comun/verificar.py` | `errores_relativos`, `verificar_fuerzas` |
| Error relativo de energía | `comun/verificar.py` | `verificar_energia` |
| Tolerancias de aprobación | `comun/verificar.py` | `TOLERANCIA_FUERZAS`, `TOLERANCIA_ENERGIA` |
| Generar todos los archivos de datos | `scripts/generar_datos.sh` | — |
