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

## Barnes-Hut en C (`c_openmp/`)

| Concepto | Archivo | Función |
|---|---|---|
| **Criterio de apertura θ** | `c_openmp/fuerzas.c` | `celda_suficientemente_lejana` |
| Interacción gravitacional con suavizado ε | `c_openmp/fuerzas.c` | `sumar_atraccion` |
| Recorrido del árbol para una partícula | `c_openmp/fuerzas.c` | `recorrer_arbol` |
| Ciclo de fuerzas sobre todas las partículas | `c_openmp/fuerzas.c` | `calcular_aceleraciones` |
| Energía cinética / potencial / total | `c_openmp/fuerzas.c` | `energia_cinetica`, `energia_potencial`, `energia_total` |
| Estructura del nodo y del árbol | `c_openmp/arbol.h` | `Nodo`, `Arbol` |
| Cubo raíz | `c_openmp/arbol.c` | `calcular_cubo_raiz` |
| Elegir el octante de una posición | `c_openmp/arbol.c` | `octante` |
| Crear los 8 hijos | `c_openmp/arbol.c` | `subdividir`, `nuevo_nodo` |
| Insertar una partícula | `c_openmp/arbol.c` | `insertar` |
| Masa y centro de masa de cada nodo | `c_openmp/arbol.c` | `calcular_masas` |
| Construcción completa (todo lo anterior) | `c_openmp/arbol.c` | `construir_arbol` |
| Partículas en el mismo punto | `c_openmp/arbol.h`, `arbol.c` | `PROFUNDIDAD_MAXIMA`, `siguiente`, `insertar` |
| Paso leapfrog completo | `c_openmp/barnes_hut.c` | `avanzar_un_paso`, `media_patada`, `deriva` |
| Cálculo inicial de fuerzas y ciclo de pasos | `c_openmp/barnes_hut.c` | `main` |
| Nombre de la versión en el CSV | `c_openmp/barnes_hut.c` | `VERSION` |
| Opciones de línea de comandos y valores por defecto | `c_openmp/configuracion.c` | `leer_configuracion` |
| Cronómetro | `c_openmp/medicion.c` | `segundos_actuales` |
| Tiempos por fase, resumen y CSV | `c_openmp/medicion.c` | `imprimir_paso`, `imprimir_resumen`, `agregar_fila_csv` |
| Lectura / escritura de archivos | `c_openmp/particulas.c` | `leer_particulas`, `escribir_particulas`, `escribir_aceleraciones` |
| Compilar normal / con sanitizadores | `c_openmp/Makefile` | `make`, `make sanitizado` |
