# Fase 3: C con OpenMP

## Qué se construyó

La versión de la fase 2 se paralelizó con OpenMP. Los cambios:

| Archivo | Qué cambió |
|---|---|
| `morton.c` / `morton.h` (nuevo) | cubo raíz con reducción paralela, códigos Morton, merge sort paralelo con tareas |
| `arbol.c` / `arbol.h` | construcción por Morton con tareas, en lugar de inserción; reserva de nodos por bloques de cada hilo |
| `fuerzas.c` / `fuerzas.h` | región paralela de fuerzas con tres repartos (`estatico`, `dinamico`, `costo`), tiempo por hilo y desbalance |
| `particulas.c` / `particulas.h` | identificador original y costo de cada partícula; reacomodo en paralelo; salida en el orden del archivo |
| `barnes_hut.c` | `recalcular_aceleraciones` (cubo → orden → árbol → fuerzas); patadas y deriva en paralelo |
| `configuracion.c` | opciones `--hilos` y `--reparto` |
| `medicion.c` | columna de desbalance por paso y promedio en el resumen |

## Cómo se usa

```bash
make -C c_openmp
c_openmp/bin/barnes_hut --entrada datos/grandes/dos_cumulos_200000.txt --pasos 10 --hilos 12 --reparto dinamico
```

`--hilos 1` es la referencia secuencial para medir la aceleración. Si no se da `--hilos`, se usan `OMP_NUM_THREADS` o todos los hilos del equipo.

## Cómo funciona

### El paso completo

```
media_patada → deriva → cubo raíz → orden Morton → árbol → fuerzas → media_patada
```

Cada flecha es el final de una región paralela, y toda región paralela termina con una **barrera implícita**: ningún hilo empieza la fase siguiente hasta que todos terminaron la actual. Así se cumple la "foto congelada": el árbol se arma cuando todas las partículas ya se movieron, y nadie se mueve hasta que todas las fuerzas están calculadas.

### 1. Cubo raíz en paralelo (`calcular_cubo_raiz`)

Cada hilo busca el mínimo y el máximo de su parte de las partículas. `reduction(min: ...)` y `reduction(max: ...)` hacen que OpenMP le dé a cada hilo su propia copia de las variables y las combine al final. Si los hilos escribieran directamente en `min_x` compartido, habría una condición de carrera.

### 2. Orden de Morton (`morton.c`)

**El código.** El cubo raíz se divide en 2²¹ rebanadas por eje. Cada coordenada se convierte en un número de 21 bits (`rebanada`), y los bits de x, y, z se intercalan (`intercalar_bits`):

```
código = ... z₁ y₁ x₁ z₀ y₀ x₀      (63 bits en total)
```

Los tres bits más altos dicen en qué hijo de la raíz está la partícula, los tres siguientes en qué nieto, y así hasta el nivel 21. El orden x, y, z de cada grupo coincide con la numeración de hijos de la fase 2: bit 0 = x, bit 1 = y, bit 2 = z.

**Por qué sirve ordenar por ese código.** Al ordenar las partículas por su código, las de cualquier celda del octree quedan **contiguas** en el arreglo: primero todas las del hijo 0, luego las del hijo 1, etc., y dentro de cada hijo, lo mismo con los nietos. Esto da tres ventajas:

1. **Construir el árbol es cortar rangos**, sin insertar una por una y sin candados.
2. **Localidad de caché**: partículas consecutivas están cerca en el espacio y recorren casi los mismos nodos del árbol. Por eso, con 1 solo hilo, esta versión ya es 2,7 veces más rápida que la de la fase 2.
3. **Repartir rangos contiguos es repartir regiones compactas del espacio**, que es lo que se usa en el reparto de trabajo (y lo que usará MPI).

**El ordenamiento** (`ordenar_mezclando`) es un merge sort con tareas: las dos mitades se ordenan en tareas paralelas, un `taskwait` espera a que ambas terminen y luego se mezclan. Por debajo de 8.192 elementos se usa `qsort` directo. Si dos partículas tienen el mismo código se desempata por índice, así el orden es siempre el mismo.

**Reacomodo físico** (`reordenar_particulas`): las partículas se copian a su nueva posición en arreglos auxiliares y luego se intercambian los punteros. Cada partícula guarda su `identificador` (su fila en el archivo), y al escribir resultados se usa para devolverlas al orden original, que es el que espera `verificar.py`.

### 3. Construcción paralela del árbol (`arbol.c`)

`construir_nodo(indice, primera, cantidad, nivel)` arma el subárbol de las partículas `[primera, primera + cantidad)`:

- **0 o 1 partícula** (o nivel 21, sin más bits que mirar): es una hoja. Se calcula su masa.
- **Si no**: se reservan 8 hijos. El rango se corta en 8 pedazos recorriéndolo y mirando los 3 bits del nivel (`octante_en_nivel`); como está ordenado, los pedazos son contiguos. Cada hijo se construye recursivamente: **en una tarea aparte si tiene más de `UMBRAL_TAREA` (1.024) partículas**, y en el mismo hilo si tiene menos.
- `#pragma omp taskwait` espera a los 8 hijos, y solo entonces se suman sus masas para obtener la del nodo.

Un hilo empieza en la raíz (`omp single`) y los demás toman las tareas a medida que aparecen. **Esto balancea la construcción solo**: un cúmulo denso genera muchas tareas, y las toma cualquier hilo que esté libre.

**Por qué no hay carreras en la construcción.** Cada subárbol escribe únicamente en sus propios nodos y lee únicamente su rango de partículas. Dos tareas nunca escriben el mismo nodo. Lo único compartido es el contador de nodos.

**El contador de nodos** (`reservar_hijos`). Dos hilos no pueden recibir los mismos 8 lugares. La primera versión hacía un `atomic capture` sobre `num_nodos` en cada subdivisión y la construcción **no escalaba**: 21 ms con 1 hilo y 16 ms con 12. Todos los hilos competían por la misma variable, y su línea de caché viajaba de núcleo en núcleo en cada una de las ~90.000 subdivisiones. Ahora cada hilo toma del contador un **bloque de 4.096 nodos** con un solo `atomic` y reparte de ahí sin sincronizar. Resultado: de 22 ms con 1 hilo a 8 ms con 8 hilos. La estructura `BloqueDeNodos` mide 64 bytes (una línea de caché) para que los bloques de dos hilos no compartan línea (**false sharing**).

El costo es que al final de cada bloque quedan nodos sin usar (a lo sumo 4.096 por hilo); por eso el programa informa `nodos_en_uso`.

**Si faltan nodos.** No se puede hacer `realloc` mientras otros hilos escriben en el arreglo. Si el contador pasa de la capacidad, `reservar_hijos` devuelve `SIN_HIJOS`, el subárbol queda incompleto y, al terminar la región paralela, `construir_arbol` duplica la capacidad y construye de nuevo. Pasa a lo sumo un par de veces en toda la simulación, porque la capacidad nueva se conserva.

**Nodo.** La hoja ya no guarda "su partícula" sino un rango: `primera_particula` y `num_particulas`. Si varias partículas tienen el mismo código (están a menos de 2⁻²¹ del lado del cubo raíz entre sí), quedan en la misma hoja y se interactúa con cada una de forma exacta. Ya no hace falta el arreglo `siguiente` de la fase 2.

### 4. Fuerzas en paralelo (`fuerzas.c`)

El recorrido es el mismo de la fase 2: mismo criterio, misma interacción. Ahora además cuenta cuántas atracciones sumó cada partícula (`interacciones`), que es su **costo**.

**Por qué no hacen falta candados.** Cada partícula lee el árbol y las posiciones (que nadie modifica en esta fase) y escribe solo `aceleraciones[i]`, `potenciales[i]` e `interacciones[i]`. Dos hilos con distinto `i` nunca escriben el mismo lugar.

**Los tres repartos** (`--reparto`):

| Reparto | Cómo | Problema que resuelve / que tiene |
|---|---|---|
| `estatico` | rangos contiguos con la misma cantidad de partículas | lo más simple; pero una partícula de un cúmulo cuesta el doble que una aislada, así que unos hilos terminan antes |
| `dinamico` (por defecto) | `schedule(dynamic, 64)`: cada hilo toma 64 partículas, y al terminar pide otras 64 | se adapta solo; el costo es pedir trabajo a una cola compartida unas 3.000 veces por paso |
| `costo` | rangos contiguos con la misma **suma de interacciones del paso anterior** (`limites_por_costo`) | sin cola compartida; supone que el costo cambia poco de un paso al siguiente |

Los límites del estático y del costo se calculan en `calcular_limites`, la función marcada con `REPARTO DE TRABAJO ENTRE HILOS`.

**Medir el desbalance.** Cada hilo mide cuánto tardó en su parte (`omp_get_wtime`), y `desbalance` devuelve `tiempo del hilo más lento / tiempo promedio`. 1,00 es perfecto; 1,20 significa que el hilo más lento tardó 20 % más que el promedio, y todos los demás lo esperaron en la barrera.

En el reparto dinámico se usa `nowait` en el `omp for`: cada hilo, al quedarse sin bloques, anota su tiempo sin esperar a los demás. Es seguro porque lo único que hace después es escribir `tiempos[hilo]`, que es solo suyo; la barrera real está al cerrar la región paralela.

**`omp_set_dynamic(0)`** (en `main`). Sin esto, OpenMP podría crear una región con menos hilos de los pedidos. En los repartos estático y por costo, las partículas del rango de un hilo que no existe no se calcularían nunca.

### 5. Integración y energía

`media_patada` y `deriva` son `parallel for`: cada `i` toca solo su partícula. Las energías usan `reduction(+ : suma)`.

## Resultados (WSL, Ryzen 5 5500U: 6 núcleos, 12 hilos)

### Correctitud

- Con θ = 0 y θ = 0,5, **el error contra la suma directa es exactamente el de la fase 2** en las tres distribuciones (0,48 %, 0,41 % y 0,36 % de promedio con θ = 0,5). El árbol por Morton tiene la misma forma que el de inserción.
- **Las aceleraciones son idénticas byte a byte** con 1, 3, 8 y 12 hilos y con los tres repartos. El resultado no depende del paralelismo: el árbol es el mismo y cada partícula suma sus atracciones en el mismo orden. Es también una prueba indirecta de que no hay carreras.
- Energía en 100 pasos con 6 hilos: igual que en la fase 2 (por ejemplo, 6,3e-6 en uniforme, según `verificar.py`).
- `make sanitizado` con 4 y 6 hilos, los tres repartos y el caso de partículas repetidas: sin avisos.

### Escalado (N = 200.000, θ = 0,5, reparto por costo, ms por paso, una corrida de 3 pasos)

| Hilos | uniforme: total | aceleración | dos_cumulos: total | aceleración | desbalance |
|---|---|---|---|---|---|
| 1 | 2465 | 1,00 | 5039 | 1,00 | 1,00 |
| 2 | 1531 | 1,61 | 2808 | 1,79 | 1,00 |
| 4 | 735 | 3,35 | 1717 | 2,93 | 1,01–1,04 |
| 6 | 667 | 3,69 | 1317 | 3,83 | 1,03 |
| 8 | 540 | 4,56 | 1202 | 4,19 | 1,02 |
| 12 | 511 | 4,82 | 960 | 5,25 | 1,02–1,09 |
| 16 | 494 | 4,99 | 1001 | 5,03 | 1,15–1,16 |

Son corridas únicas y cortas, así que tienen ruido. Las medianas de 3–5 corridas se hacen en la fase 7.

**Comparado con la fase 2**: el paso uniforme de 200.000 bajó de 6,67 s a 2,47 s con 1 hilo (solo por el orden de Morton) y a 0,51 s con 12 hilos (13 veces más rápido).

**Lectura de las curvas:**

- **De 1 a 6 hilos** se gana casi 4 veces: hay 6 núcleos físicos.
- **De 6 a 12 hilos** se gana mucho menos (1,3–1,4 veces): los hilos 7 a 12 son los segundos hilos SMT de cada núcleo. Comparten las unidades de cálculo y la caché con el primero, así que solo aprovechan los ratos en que el otro espera a la memoria.
- **Con 16 hilos** no se gana nada: hay más hilos que hilos de hardware, el sistema operativo los turna y el desbalance sube a 1,15 porque algunos hilos pierden su turno a mitad del trabajo.
- **De 1 a 2 hilos** se gana menos de 2 veces. Una hipótesis para comprobar en la fase 7: con un solo núcleo activo el procesador sube más su frecuencia (turbo), y con más núcleos activos la baja.

**Fases que no escalan.** Con 8 hilos, la construcción (cubo + árbol) tarda unos 9 ms y el orden unos 13 ms por paso, frente a 520 ms de fuerzas: cerca del 4 % del paso. El orden escala poco porque la última mezcla del merge sort la hace un solo hilo, y porque reacomodar 7 arreglos de 200.000 elementos está limitado por el ancho de banda de memoria, no por los núcleos. Por la ley de Amdahl, una fracción secuencial del 4 % limita la aceleración máxima a 25.

### Con N pequeño, más hilos no ayudan tanto (N = 2.000, dos_cumulos, 100 pasos)

| Fase (ms por paso) | 1 hilo | 12 hilos |
|---|---|---|
| construcción | 0,28 | 1,29 |
| ordenamiento | 0,30 | 1,96 |
| fuerzas | 25,31 | 4,40 |
| integración | 0,02 | 1,33 |
| total | 25,91 | 8,98 |

Las fuerzas se aceleran 5,8 veces, pero el orden, la construcción y la integración se vuelven **más lentas** con 12 hilos. Cada región paralela cuesta despertar y sincronizar a los hilos (microsegundos), y con 2.000 partículas el trabajo de esas fases es menor que ese costo. Pasan de ser el 2 % del paso a ser el 51 %, y la aceleración total queda en 2,9. El desbalance de las fuerzas también sube (1,26), porque con tan pocas partículas los bloques de 64 son pocos y un hilo que tarda en arrancar ya se nota.

### Comparación de repartos (N = 200.000, tiempo de fuerzas por paso)

| Distribución | Hilos | estático | dinámico | costo |
|---|---|---|---|---|
| dos_cumulos | 6 | 1423 ms (desb. 1,05) | 1358 ms (1,00) | 1392 ms (1,02) |
| dos_cumulos | 12 | 1026 ms (1,10) | 923 ms (1,00) | 936 ms (1,02) |
| concentrada | 6 | 1481 ms (1,18) | 1301 ms (1,00) | 1367 ms (1,03) |
| concentrada | 12 | 974 ms (1,20) | 850 ms (1,00) | 904 ms (1,03) |

- **El estático pierde más con `concentrada`** (desbalance 1,20): los hilos a los que les toca el cúmulo tienen partículas más caras.
- **El dinámico es el más rápido** y por eso es el reparto por defecto.
- **El reparto por costo queda cerca** (2–6 % más lento). Su predicción usa el paso anterior y no ve los efectos de caché. Se mantiene porque es la idea que se necesita en MPI, donde no hay una cola compartida de la que los procesos puedan pedir trabajo.

### Qué se ve en htop

Durante las fuerzas, las 12 barras (o las que se pidan con `--hilos`) en 100 %. Entre paso y paso hay un instante en que bajan, durante el orden y la construcción. Con `--hilos 16`, htop muestra el programa usando todos los hilos de hardware, pero el tiempo por paso no mejora respecto a 12.

## Dónde se hacen las modificaciones en vivo

| Pedido | Dónde | Qué cambiar |
|---|---|---|
| Cambiar el reparto | línea de comandos | `--reparto estatico`, `dinamico` o `costo` |
| Cambiar cómo se reparte (por ejemplo, intercalado: el hilo h toma h, h+H, h+2H…) | `fuerzas.c`, `calcular_aceleraciones` | agregar una rama con `schedule(static, 1)` |
| Otro criterio de costo (por ejemplo, tiempo medido en vez de interacciones) | `fuerzas.c`, `limites_por_costo` | usar otro arreglo como peso |
| Tamaño de bloque del dinámico | `fuerzas.c`, `BLOQUE_DINAMICO` | — |
| Criterio θ | `fuerzas.c`, `celda_suficientemente_lejana` | igual que en la fase 2 |
| Granularidad de las tareas del árbol | `arbol.c`, `UMBRAL_TAREA` | más grande: menos tareas y menos balanceo |

## Decisiones de diseño

- **Morton y no candados.** Con candados (uno por nodo), cada inserción bloquea los nodos que modifica. Funciona, pero cerca de la raíz todos los hilos compiten por los mismos candados, el orden de inserción cambia de una corrida a otra (y con él el árbol) y hay que razonar con cuidado para no caer en interbloqueos. Con Morton no hay candados, el árbol siempre es el mismo y, además, se gana la localidad de caché.
- **Tareas y no `parallel for` para el árbol.** El árbol es irregular: un cúmulo tiene subárboles mucho más grandes que el vacío. Las tareas reparten el trabajo a medida que aparece.
- **Bloques de nodos por hilo.** El `atomic` por subdivisión era correcto, pero no escalaba. Es un ejemplo de que una operación atómica no es gratis: serializa a los hilos que tocan la misma variable.
- **Reacomodo físico de las partículas**, no solo un arreglo de índices: así el ciclo de fuerzas lee posiciones y masas consecutivas en memoria.
- **El reparto dinámico por defecto**, porque fue el más rápido medido. El estático queda para mostrar el problema y el de costo para preparar MPI.

## Preguntas tipo sustentación

1. **¿Por qué el cálculo de fuerzas no necesita candados ni `atomic`?**
   Cada partícula lee datos que nadie modifica en esa fase (árbol y posiciones) y escribe solo en su propia posición `i` de los arreglos de resultados. No hay dos hilos que escriban el mismo lugar.

2. **¿Qué pasaría si se quita el `taskwait` de `construir_nodo`?**
   El nodo sumaría las masas de sus hijos antes de que las tareas de los hijos terminen: leería masas en 0 o a medio calcular. Es una condición de carrera. El resultado cambiaría de una corrida a otra y las fuerzas serían incorrectas.

3. **¿Por qué la primera versión de la construcción no escalaba y cómo se arregló?**
   Cada subdivisión hacía un `atomic` sobre el mismo contador. El `atomic` es correcto, pero todos los hilos pelean por la misma línea de caché y se serializan. Se arregló reservando bloques de 4.096 nodos por hilo: un `atomic` por bloque en lugar de uno por subdivisión.

4. **¿Qué es false sharing y dónde se evitó?**
   Dos hilos escriben variables distintas que están en la misma línea de caché (64 bytes). El hardware trata la línea como compartida y la pasa de un núcleo a otro en cada escritura. Se evitó haciendo que `BloqueDeNodos` mida 64 bytes.

5. **¿Qué pasa al pasar de 6 a 12 hilos? ¿Y a 16?**
   De 6 a 12 se gana poco (1,3–1,4 veces), porque los hilos 7 a 12 son SMT y comparten núcleo con los primeros 6. Con 16 no se gana nada, porque hay más hilos que hilos de hardware: el sistema operativo los turna y el desbalance sube.

6. **¿Por qué el reparto estático es el peor con la distribución concentrada?**
   Reparte la misma cantidad de partículas, pero no el mismo trabajo. Las partículas del cúmulo abren más celdas y cuestan más, así que los hilos que reciben el cúmulo tardan hasta 20 % más, y el resto los espera en la barrera.

7. **¿Cómo funciona el reparto por costo y en qué supuesto se basa?**
   Cuenta las interacciones de cada partícula en el paso anterior, y corta el arreglo (ordenado por Morton) en rangos contiguos con la misma suma. Supone que las partículas se mueven poco entre pasos, así que el costo del paso anterior predice bien el actual.

8. **¿Por qué ordenar por Morton acelera incluso con un solo hilo?**
   Porque partículas consecutivas en el arreglo están cerca en el espacio y recorren casi los mismos nodos del árbol, que entonces ya están en la caché. Con el orden aleatorio del archivo, casi cada nodo había que traerlo de la RAM.

9. **¿Por qué se usa `omp_set_dynamic(0)`?**
   Para que cada región paralela tenga exactamente el número de hilos pedido. En los repartos estático y por costo, si OpenMP creara menos hilos, los rangos de los hilos que faltan no se calcularían y esas partículas quedarían con la aceleración del paso anterior.

10. **Si las fuerzas se aceleran 8 veces y el orden y la construcción no, ¿qué limita la aceleración total?**
    La ley de Amdahl: la parte que no se acelera pesa cada vez más. Con 8 hilos, orden y construcción son cerca del 4 % del paso; con infinitos hilos para las fuerzas, la aceleración nunca pasaría de 1 / 0,04 = 25.

11. **¿Por qué las aceleraciones dan idénticas byte a byte con cualquier número de hilos? ¿Pasa lo mismo con la energía?**
    Las aceleraciones de cada partícula se suman siempre en el mismo orden, porque el árbol es el mismo. Por eso, después de 100 pasos, el estado final con 1 y con 12 hilos es idéntico. La energía es una reducción: cada hilo suma su parte y luego se combinan. Como la suma en coma flotante no es asociativa, otro agrupamiento puede cambiar los últimos bits. En la práctica coincide en los 10 decimales que se imprimen.

12. **Con 2.000 partículas, 12 hilos dan solo 2,9 veces de aceleración. ¿Por qué?**
    Porque hay poco trabajo. Cada región paralela tiene un costo fijo para despertar y sincronizar a los hilos. Con N pequeño, en el orden, la construcción y la integración ese costo es mayor que el trabajo, y esas fases se vuelven más lentas que con 1 hilo. Es la misma ley de Amdahl, con una parte "secuencial" que crece con el número de hilos.
