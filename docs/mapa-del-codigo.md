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

## C con OpenMP (`c_openmp/`)

### Lo que se puede pedir modificar en vivo

| Concepto | Archivo | Función / constante |
|---|---|---|
| **Criterio de apertura θ** | `c_openmp/fuerzas.c` | `celda_suficientemente_lejana` |
| **Reparto de partículas entre hilos** (estático y por costo) | `c_openmp/fuerzas.c` | `calcular_limites`, `limites_por_costo` |
| Reparto dinámico (tamaño de bloque) | `c_openmp/fuerzas.c` | `BLOQUE_DINAMICO`, `calcular_aceleraciones` |
| Reparto por defecto | `c_openmp/configuracion.c` | `.reparto` en `leer_configuracion` |
| Tamaño mínimo de un subárbol para crear una tarea | `c_openmp/arbol.c` | `UMBRAL_TAREA` |
| Tamaño mínimo para ordenar en paralelo | `c_openmp/morton.c` | `UMBRAL_ORDENAMIENTO_SECUENCIAL` |

### Fases del paso

| Concepto | Archivo | Función |
|---|---|---|
| Paso leapfrog completo | `c_openmp/barnes_hut.c` | `avanzar_un_paso`, `media_patada`, `deriva` |
| Cubo raíz → orden → árbol → fuerzas (con tiempos) | `c_openmp/barnes_hut.c` | `recalcular_aceleraciones` |
| Cálculo inicial de fuerzas y ciclo de pasos | `c_openmp/barnes_hut.c` | `main` |
| Número de hilos fijo (`omp_set_dynamic(0)`) | `c_openmp/barnes_hut.c` | `main` |

### Ordenamiento por Morton

| Concepto | Archivo | Función |
|---|---|---|
| Cubo raíz (reducción paralela de mínimos y máximos) | `c_openmp/morton.c` | `calcular_cubo_raiz` |
| Coordenada → rebanada de 21 bits | `c_openmp/morton.c` | `rebanada` |
| Intercalar bits x, y, z (código Morton) | `c_openmp/morton.c` | `intercalar_bits`, `calcular_codigos` |
| Merge sort paralelo con tareas | `c_openmp/morton.c` | `ordenar_mezclando`, `mezclar`, `comparar_pares` |
| Ordenar y reacomodar las partículas | `c_openmp/morton.c` | `ordenar_por_morton` |
| Mover las partículas a su nuevo lugar | `c_openmp/particulas.c` | `reordenar_particulas` |
| Bits por eje | `c_openmp/morton.h` | `BITS_POR_EJE` |

### Árbol

| Concepto | Archivo | Función |
|---|---|---|
| Estructura del nodo y del árbol | `c_openmp/arbol.h` | `Nodo`, `Arbol`, `BloqueDeNodos` |
| Construcción paralela (tareas) | `c_openmp/arbol.c` | `construir_arbol`, `construir_nodo` |
| Hijo en el que cae un código en cada nivel | `c_openmp/arbol.c` | `octante_en_nivel` |
| Reserva de nodos por bloques de cada hilo (atómica) | `c_openmp/arbol.c` | `reservar_hijos`, `NODOS_POR_BLOQUE` |
| Masa y centro de masa | `c_openmp/arbol.c` | `calcular_masa_de_hoja`, `sumar_masas_de_hijos` |
| Ampliar el arreglo si faltan nodos | `c_openmp/arbol.c` | final de `construir_arbol` |

### Fuerzas y energía

| Concepto | Archivo | Función |
|---|---|---|
| Interacción gravitacional con suavizado ε | `c_openmp/fuerzas.c` | `sumar_atraccion` |
| Recorrido del árbol para una partícula | `c_openmp/fuerzas.c` | `recorrer_arbol`, `fuerza_sobre_particula` |
| Región paralela de fuerzas y tiempo por hilo | `c_openmp/fuerzas.c` | `calcular_aceleraciones` |
| Desbalance entre hilos | `c_openmp/fuerzas.c` | `desbalance` |
| Energía cinética / potencial / total | `c_openmp/fuerzas.c` | `energia_cinetica`, `energia_potencial`, `energia_total` |

### Entrada, salida y medición

| Concepto | Archivo | Función |
|---|---|---|
| Opciones de línea de comandos y valores por defecto | `c_openmp/configuracion.c` | `leer_configuracion` |
| Lectura / escritura (en el orden del archivo) | `c_openmp/particulas.c` | `leer_particulas`, `escribir_particulas`, `escribir_aceleraciones`, `posiciones_por_fila` |
| Cronómetro | `c_openmp/medicion.c` | `segundos_actuales` |
| Tiempos por fase, resumen y CSV | `c_openmp/medicion.c` | `imprimir_paso`, `imprimir_resumen`, `agregar_fila_csv` |
| Compilar normal / con sanitizadores | `c_openmp/Makefile` | `make`, `make sanitizado` |

## Python con multiprocessing (`python_mp/`)

### Lo que se puede pedir modificar en vivo

| Concepto | Archivo | Función / constante |
|---|---|---|
| **Criterio de apertura θ** | `python_mp/fuerzas.py` | `celdas_suficientemente_lejanas` |
| **Reparto de grupos entre procesos** | `python_mp/trabajadores.py` | `repartir_grupos` |
| Grupos por tarea en el reparto dinámico | `python_mp/trabajadores.py` | `GRUPOS_POR_TAREA` |
| Cómo llega el árbol a los trabajadores | `python_mp/trabajadores.py` | `obtener_arbol` (y `--estrategia`) |
| Tamaño de grupo | línea de comandos | `--tamano-grupo` (16) |
| Valores por defecto de las opciones | `python_mp/barnes_hut.py` | `leer_configuracion` |

### Árbol (arreglos planos de numpy)

| Concepto | Archivo | Función |
|---|---|---|
| Campos del árbol (un arreglo por dato del nodo) | `python_mp/arbol.py` | `Arbol` |
| Cubo raíz | `python_mp/arbol.py` | `calcular_cubo_raiz` |
| Códigos Morton | `python_mp/arbol.py` | `codigos_morton` |
| Construcción nivel por nivel | `python_mp/arbol.py` | `construir_arbol` |
| Ordenar y reacomodar las partículas | `python_mp/barnes_hut.py` | `Simulacion.ordenar_por_morton` |

### Fuerzas

| Concepto | Archivo | Función |
|---|---|---|
| Rangos de índices a un solo arreglo | `python_mp/fuerzas.py` | `expandir_rangos` |
| Recorrido del árbol para un grupo (lista de interacción) | `python_mp/fuerzas.py` | `lista_de_interaccion` |
| Matriz grupo × fuentes, sin autointeracción | `python_mp/fuerzas.py` | `fuerzas_de_grupo` |

### Procesos y memoria compartida

| Concepto | Archivo | Función |
|---|---|---|
| Arreglos de numpy en memoria compartida | `python_mp/memoria_compartida.py` | `ArreglosCompartidos.crear`, `conectar`, `cerrar` |
| Crear el pool (forkserver) y la memoria compartida | `python_mp/trabajadores.py` | `EquipoDeTrabajadores.__init__` |
| Copiar el árbol a memoria compartida | `python_mp/trabajadores.py` | `EquipoDeTrabajadores.publicar_arbol` |
| Repartir, lanzar tareas y medir desbalance | `python_mp/trabajadores.py` | `EquipoDeTrabajadores.calcular_aceleraciones` |
| Inicio de cada trabajador | `python_mp/trabajadores.py` | `iniciar_trabajador` |
| Lo que hace un trabajador con una tarea | `python_mp/trabajadores.py` | `calcular_tarea` |
| Liberar procesos y memoria al final | `python_mp/trabajadores.py`, `barnes_hut.py` | `EquipoDeTrabajadores.cerrar`, `main` |

### Paso, medición y archivos

| Concepto | Archivo | Función |
|---|---|---|
| Hilos de numpy fijados en 1 | `python_mp/barnes_hut.py` | inicio del archivo |
| Paso leapfrog completo | `python_mp/barnes_hut.py` | `avanzar_un_paso`, `media_patada`, `deriva` |
| Orden → árbol → fuerzas (con tiempos) | `python_mp/barnes_hut.py` | `recalcular_aceleraciones` |
| Energía | `python_mp/barnes_hut.py` | `Simulacion.energia_total` |
| Resultados en el orden del archivo | `python_mp/barnes_hut.py` | `Simulacion.en_orden_del_archivo` |
| Tiempos, resumen y CSV | `python_mp/medicion.py` | `TiemposPorFase`, `imprimir_resumen`, `agregar_fila_csv` |
