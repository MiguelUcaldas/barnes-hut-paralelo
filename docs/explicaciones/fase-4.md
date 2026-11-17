# Fase 4: Python con multiprocessing

## Qué se construyó

Una versión completa en `python_mp/`, con numpy y `multiprocessing` y sin otras librerías.

| Archivo | Contenido |
|---|---|
| `barnes_hut.py` | opciones, estado de las partículas (`Simulacion`), orden Morton, leapfrog, `main` |
| `arbol.py` | `Arbol` (un arreglo de numpy por cada dato del nodo), cubo raíz, códigos Morton, construcción nivel por nivel |
| `fuerzas.py` | **criterio de apertura**, lista de interacción de un grupo, matriz grupo × fuentes |
| `trabajadores.py` | pool de procesos, **las tres estrategias** para que el árbol llegue a los trabajadores, **reparto de grupos** |
| `memoria_compartida.py` | `ArreglosCompartidos`: arreglos de numpy sobre bloques de memoria compartida |
| `medicion.py` | tiempos por fase, resumen y fila del CSV (mismo formato que la versión en C) |

## Cómo se usa

```bash
.venv/bin/python python_mp/barnes_hut.py --entrada datos/grandes/dos_cumulos_20000.txt --pasos 10 \
    --trabajadores 6 --estrategia compartida --reparto dinamico --tamano-grupo 16
```

Valores por defecto: `--estrategia compartida`, `--reparto dinamico`, `--tamano-grupo 16`, todos los núcleos.

## El problema de Python y cómo se resolvió

Un ciclo de Python ejecuta del orden de 10 millones de operaciones simples por segundo. Barnes-Hut con 200.000 partículas hace del orden de 10⁹ interacciones por paso. Recorrer el árbol partícula por partícula en Python tardaría horas por paso.

La solución es que **Python nunca haga un ciclo por partícula ni por nodo**: los ciclos de Python recorren niveles del árbol o grupos de partículas, y todo el trabajo por partícula o por nodo lo hace numpy, en C, sobre arreglos completos.

### 1. El árbol en arreglos planos (`arbol.py`)

No hay un objeto por nodo. Un nodo es un índice `k`, y sus datos están en `tamano[k]`, `masa[k]`, `centro_de_masa[k]`, `primer_hijo[k]`, `num_hijos[k]`, `primera_particula[k]` y `num_particulas[k]`.

A diferencia de la versión en C, **solo existen los nodos con partículas**: un nodo tiene de 1 a 8 hijos, guardados seguidos desde `primer_hijo`.

**Construcción nivel por nivel.** Las partículas se ordenan por código Morton (igual que en la fase 3). En el nivel L, dos partículas están en la misma celda si sus códigos comparten los primeros 3·L bits; como están ordenadas, cada celda es un **tramo contiguo** del arreglo. Para cada nivel:

1. `prefijos = codigos >> (3·(21 − L))` para todas las partículas a la vez.
2. Las celdas empiezan donde el prefijo cambia: `prefijos[1:] != prefijos[:-1]`.
3. Una celda es un nodo nuevo solo si su padre tenía 2 o más partículas. Si tenía una, el padre ya era hoja.

Son a lo sumo 21 iteraciones de Python, cada una con unas diez operaciones de numpy sobre todo el arreglo.

Los nodos quedan ordenados por nivel y, dentro de cada nivel, por padre, así que los hijos de un mismo padre quedan seguidos. `np.unique` encuentra el primer hijo y la cantidad de hijos de cada padre.

**Masa y centro de masa sin ciclos.** Con las sumas acumuladas `M[i] = m₀ + … + m_{i−1}`, la masa del tramo `[a, b)` es `M[b] − M[a]`: una resta por nodo, todas a la vez. Igual con el momento Σ m·x, que dividido por la masa da el centro de masa.

### 2. Fuerzas por grupos (`fuerzas.py`)

Se toman grupos de 16 partículas consecutivas en orden Morton (están cerca en el espacio), y para cada grupo:

1. **Lista de interacción** (`lista_de_interaccion`): se recorre el árbol nivel por nivel. La `frontera` tiene todos los nodos de un nivel que hay que revisar, y el criterio se evalúa para todos a la vez:
   - las celdas lejanas van a la lista;
   - las hojas cercanas aportan sus partículas;
   - los nodos internos cercanos se abren, y sus hijos forman la frontera siguiente.
2. **Matriz grupo × fuentes** (`fuerzas_de_grupo`): con las 16 partículas y las K fuentes de la lista (celdas y partículas) se arman matrices de 16 × K y se calcula todo de una vez. `dx`, `dy`, `dz` → `1/√(r² + ε²)` → aceleración y potencial.

**El criterio por grupo** (`celdas_suficientemente_lejanas`) es el mismo `tamaño / distancia < θ`, pero la distancia se mide desde el centro de masa de la celda **hasta el punto más cercano de la caja del grupo**. Si la celda está lejos de ese punto, está lejos de todas las partículas del grupo. Por eso es un poco más estricto que el criterio por partícula: da menos error con el mismo θ, a cambio de más interacciones. Con θ = 0 sigue siendo la suma directa.

**Autointeracción.** Las hojas de las partículas del grupo nunca están lejos del grupo (la distancia a la caja es 0), así que esas partículas aparecen en su propia lista. Su posición en la matriz se pone en cero.

`expandir_rangos` convierte rangos `[inicio, inicio + cantidad)` en un solo arreglo de índices sin ciclos. Se usa para pasar de "abrir estos nodos" a "estos son sus hijos" y de "estas hojas" a "estas partículas".

### 3. Procesos (`trabajadores.py`)

Al iniciar se crea un `Pool` de procesos con el método **`forkserver`**. Cada trabajador arranca desde un proceso limpio y no hereda la memoria del principal: todo lo que comparten es explícito. Con `fork` los hijos verían una copia de los arreglos del padre "gratis", y la comparación entre estrategias no sería honesta.

En cada paso:

1. El **principal** hace la patada y la deriva, calcula el cubo raíz y ordena por Morton (en memoria compartida).
2. Según la estrategia, construye el árbol y lo publica.
3. Reparte los grupos en tareas y llama a `pool.map`, que **no regresa hasta que todas las tareas terminaron**. Esa es la barrera de la "foto congelada": nadie mueve partículas mientras haya fuerzas por calcular.
4. Cada **trabajador** calcula sus grupos y escribe las aceleraciones directamente en memoria compartida. Cada grupo es un tramo distinto del arreglo, así que no hay dos procesos escribiendo en el mismo lugar.

**Memoria compartida** (`memoria_compartida.py`). `SharedMemory` crea un bloque de memoria del sistema operativo con un nombre. `np.ndarray(forma, buffer=bloque.buf)` lo usa como arreglo de numpy sin copiar. A los trabajadores solo se les envía la descripción (nombre, forma y tipo), que son unos bytes. Al final, el principal borra los bloques (`unlink`). Por eso `main` usa `try/finally`: si no se borran, sobreviven al programa en `/dev/shm`.

**Hilos de numpy en 1.** Al inicio de `barnes_hut.py` se fija `OMP_NUM_THREADS=1` (y las variables equivalentes) antes de importar numpy. Si numpy usara varios hilos internos, 12 procesos × varios hilos serían muchos más hilos que núcleos.

## La decisión que pide el enunciado: ¿cómo llega el árbol a los procesos?

Se implementaron las tres opciones del enunciado (`--estrategia`). En las tres, las partículas y los resultados están en memoria compartida; solo cambia el árbol:

| Estrategia | Cómo funciona |
|---|---|
| `compartida` | el principal construye el árbol una vez y lo copia a memoria compartida; los trabajadores lo leen sin copiarlo |
| `local` | cada trabajador construye su propio árbol a partir de las partículas compartidas (trabajo repetido, cero comunicación) |
| `copia` | el principal construye el árbol y lo envía serializado con pickle dentro de cada tarea (comunicación explícita) |

### Mediciones (N = 200.000; el árbol tiene unos 295.000 nodos y pesa 21 MB)

**Costo de hacer llegar el árbol, medido aparte:**

| | 1 trabajador | 6 trabajadores | 12 trabajadores |
|---|---|---|---|
| `compartida`: construir una vez + copiar a memoria compartida | 160 + 4 ms | 160 + 4 ms | 160 + 4 ms |
| `local`: la construcción más lenta entre P construcciones simultáneas | 130 ms | 190–220 ms | 290–380 ms |
| `copia`: construir + enviar P copias por las tuberías del Pool | 160 + 75 ms | 160 + 330 ms | 160 + 490–565 ms |
| memoria usada por los árboles | 21 MB | 21 MB / 126 MB / 126 MB | 21 MB / 252 MB / 252 MB |

(En la última fila, el orden es compartida / local / copia.)

**Tiempo total por paso** (uniforme, 2 pasos): entre 4,15 y 5,0 s con cualquier estrategia, con diferencias de 2–10 % que están dentro del ruido. Las fuerzas son el 95 % del paso y tapan la diferencia. En `copia` la transferencia también se solapa en parte con el cálculo: el primer trabajador empieza mientras se envía la copia al segundo.

**Decisión: memoria compartida.**

- Es la única estrategia cuyo costo **no crece con el número de trabajadores**: un solo árbol, construido una vez y copiado en 4 ms.
- En `local`, las P construcciones simultáneas compiten por el ancho de banda de memoria y tardan más cuantos más trabajadores hay (de 130 a ~380 ms), además de ocupar P veces la memoria.
- En `copia`, el principal tiene que serializar y enviar P copias de 21 MB, uno tras otro. Crece linealmente con P, también ocupa P veces la memoria, y con reparto dinámico serían cientos de copias por paso (por eso el programa no lo permite).
- Con N grande el árbol pesa más y la diferencia crece. Con MPI (fase 5) no hay memoria compartida y esta opción no existe: por eso allá el enunciado pide enviar solo las partes del árbol que cada proceso necesita.

## Resultados

### Correctitud

| θ | uniforme (prom. / máx.) | dos_cumulos | concentrada |
|---|---|---|---|
| 0 | 1,4e-10 / 4,2e-10 | 7,3e-11 / 4,8e-10 | 9,1e-11 / 4,6e-10 |
| 0,5 | 0,16 % / 5,0 % | 0,090 % / 1,5 % | 0,081 % / 2,5 % |

Con grupos de 16 y θ = 0,5. El error es 3 a 5 veces menor que en C (0,36–0,48 %) por el criterio por grupo, que es más estricto. Con grupos de 64 es aún menor (0,087 % en uniforme), porque la caja del grupo es más grande y se abren más celdas.

- **Resultados idénticos byte a byte** con las tres estrategias, los tres repartos y 1, 3, 5, 8 o 12 trabajadores. El árbol y los grupos no dependen de cuántos procesos haya.
- **Energía** en 100 pasos con 6 trabajadores, según `verificar.py`: 5,6e-6 (uniforme), 3,2e-5 (dos_cumulos), 2,0e-5 (concentrada). Todas pasan.
- **Casos límite** (partículas repetidas, todas en un punto, una sola partícula): exactos con θ = 0.
- **Memoria compartida**: después de cada corrida `/dev/shm` queda vacío.

### Tamaño de grupo (N = 20.000, dos_cumulos, tiempo de fuerzas por paso)

| Grupo | 1 trabajador | 6 trabajadores | 12 trabajadores |
|---|---|---|---|
| 8 | 3036 ms | 857 ms | 754 ms |
| 16 | 1973 ms | 649 ms | 559 ms |
| 24 | 1602 ms | 631 ms | 676 ms |
| 32 | 1642 ms | 866 ms | 826 ms |
| 64 | 1999 ms | 1307 ms | — |
| 256 | 3370 ms | 2544 ms | — |

- **Con grupos chicos** pesa el costo de Python: hay muchos grupos y cada recorrido tiene unas 15 llamadas a numpy por nivel.
- **Con grupos grandes**, cada matriz de grupo × fuentes ocupa varios MB. No cabe en la caché, y cada operación de numpy lee y escribe la RAM. Con un trabajador da igual, pero con 6 todos compiten por el mismo ancho de banda de memoria: con grupos de 64 se gana solo 1,5 veces, con grupos de 16 se gana 3 veces.
- Se probó también partir la lista de fuentes en bloques de 1.024 o 4.096 para que cupieran en la caché L2, y fue más lento (más llamadas a numpy). Se descartó.

### Escalado (N = 20.000, dos_cumulos, compartida, dinámico, ms por paso)

| Trabajadores | 1 | 2 | 4 | 6 | 8 | 12 | 16 |
|---|---|---|---|---|---|---|---|
| total | 2004 | 1175 | 698 | 588 | 603 | 603 | 583 |
| aceleración | 1,00 | 1,71 | 2,87 | 3,41 | 3,32 | 3,32 | 3,44 |

Con N = 200.000 (uniforme): 14,5 s por paso con 1 trabajador y ~4,2 s con 6 o 12 (3,4 veces).

**Por qué la curva se aplana en 6** (en C seguía subiendo hasta 12):

- Las operaciones de numpy recorren arreglos temporales enteros (`dx`, `dy`, `dz`, `inversa`, `peso`…), y cada operación lee y escribe memoria. El trabajo está más limitado por el **ancho de banda de memoria** que por el cálculo.
- Los dos hilos SMT de un núcleo comparten las cachés L1 y L2 y el mismo acceso a memoria, así que el segundo hilo casi no tiene con qué ayudar.
- La versión en C, en cambio, hace las cuentas en registros y depende más del cálculo.
- Con 16 trabajadores tampoco empeora, porque el cuello de botella sigue siendo la memoria y no la CPU.

**Fracción secuencial.** El orden y la construcción en el principal tardan unos 15 ms con N = 20.000 (2,5 % del paso con 6 trabajadores) y unos 180 ms con N = 200.000 (4 %). Mientras tanto, los trabajadores esperan.

### Repartos (N = 20.000, concentrada, compartida)

| Trabajadores | estático | dinámico | costo |
|---|---|---|---|
| 6 | 580 ms (desb. 1,12) | 603 ms (1,01) | 591 ms (1,08) |
| 12 | 592 ms (1,17) | 574 ms (1,02) | 578 ms (1,08) |

- El estático es el más desbalanceado (1,17 con 12), como en C.
- Las diferencias de tiempo son pequeñas porque el límite es la memoria: aunque un trabajador espere, los demás no van más rápido por tener la memoria para ellos.
- El dinámico queda por defecto porque tiene el menor desbalance y fue el más rápido con 12 trabajadores.

### Python frente a C (N = 200.000 uniforme, por paso)

| | 1 hilo / trabajador | 6 | 12 |
|---|---|---|---|
| C con OpenMP | 2,47 s | 0,67 s | 0,51 s |
| Python | 14,5 s | ~4,3 s | ~4,1 s |

Python es unas 6 veces más lento con un trabajador, y 8 veces más con 12, porque escala menos.

## Qué se ve en htop

Un proceso por trabajador (todos `python`), en 100 % durante las fuerzas. Antes de cada paso, unos milisegundos con solo el principal activo (orden y construcción). La memoria compartida aparece en la columna SHR de cada proceso.

Con `--estrategia local`, al empezar cada paso todos los trabajadores construyen su árbol a la vez, y la memoria total (RES) es mayor.

## Dónde se hacen las modificaciones en vivo

| Pedido | Dónde |
|---|---|
| Criterio θ | `fuerzas.py`, `celdas_suficientemente_lejanas`. Por ejemplo, para usar el criterio por partícula como en C, la distancia tendría que medirse a cada partícula y no a la caja (la matriz de criterio pasaría a ser grupo × nodos) |
| Reparto de trabajo | `trabajadores.py`, `repartir_grupos` (o `--reparto`) |
| Tamaño de las tareas dinámicas | `trabajadores.py`, `GRUPOS_POR_TAREA` |
| Cómo llega el árbol | `--estrategia`; el código está en `obtener_arbol` y `publicar_arbol` |
| Tamaño de grupo | `--tamano-grupo` |

## Decisiones de diseño

- **Grupos y listas de interacción** en lugar de recorrer partícula por partícula: es la única forma de que numpy haga el trabajo. El precio es un criterio un poco más estricto (más interacciones, menos error).
- **El árbol solo con nodos no vacíos.** En numpy no hay costo extra por tener de 1 a 8 hijos, y el árbol pesa menos de la mitad que con 8 hijos fijos (295.000 nodos en lugar de 750.000).
- **Memoria compartida para el árbol**, por las mediciones de arriba.
- **`forkserver`**, para que lo compartido sea explícito y la comparación entre estrategias sea justa.
- **Grupos de 16**, por la medición del tamaño de grupo: es el mejor con varios trabajadores.

## Preguntas tipo sustentación

1. **¿Por qué no se recorre el árbol partícula por partícula, como en C?**
   Porque un ciclo de Python es unas 100 veces más lento que uno de C. Con 10⁹ interacciones por paso serían horas. Recorriendo por grupos, Python hace un ciclo por nivel y por grupo, y numpy hace el trabajo por partícula y por nodo.

2. **¿Cómo se construye el árbol sin un ciclo por nodo?**
   Con las partículas ordenadas por Morton, cada celda de cada nivel es un tramo contiguo. Para un nivel se calculan los prefijos de todos los códigos, se buscan los lugares donde cambian y se descartan las celdas cuyo padre tenía una sola partícula. Son 21 niveles como máximo, cada uno con operaciones de numpy sobre todo el arreglo.

3. **¿Qué estrategia eligieron para compartir el árbol y por qué?**
   Memoria compartida. Su costo (construir una vez y copiar en 4 ms) no crece con el número de trabajadores. La reconstrucción local repite la construcción en cada proceso y, como compiten por la memoria, tarda más cuantos más son (de 130 a ~380 ms). Enviar copias crece linealmente (de 75 a 565 ms). Las dos ocupan P veces la memoria.

4. **Si el árbol está en memoria compartida y todos lo leen, ¿hay condiciones de carrera?**
   No, porque durante las fuerzas nadie lo modifica: el principal lo escribe antes de `pool.map` y los trabajadores solo lo leen. Las aceleraciones se escriben en memoria compartida, pero cada trabajador escribe solo los tramos de sus grupos, que no se cruzan.

5. **¿Dónde está la barrera entre fases?**
   En `pool.map`: no regresa hasta que todas las tareas terminaron. Recién entonces el principal hace la segunda media patada y el paso siguiente.

6. **¿Por qué el error con θ = 0,5 es menor que en C?**
   Porque el criterio se evalúa contra el punto del grupo más cercano a la celda. Una celda que pasaría el criterio para una partícula del grupo puede no pasarlo para la caja, y entonces se abre. Hay más interacciones y menos error.

7. **¿Por qué la versión de Python no mejora de 6 a 12 trabajadores y la de C sí?**
   numpy trabaja con arreglos temporales enteros: cada operación lee y escribe memoria, así que el límite es el ancho de banda de memoria. Los hilos SMT comparten núcleo y caché, y no agregan ancho de banda. El código C hace las cuentas en registros y saca más provecho del segundo hilo.

8. **¿Qué pasa si el programa termina con un error y no borra la memoria compartida?**
   Los bloques quedan en `/dev/shm` ocupando RAM hasta que alguien los borre o se reinicie el sistema. Por eso `main` cierra todo en un `finally`.

9. **¿Por qué se fija `OMP_NUM_THREADS=1` antes de importar numpy?**
   Para que cada proceso use un solo hilo. Si numpy abriera sus propios hilos, con 12 procesos habría muchos más hilos que núcleos, y el escalado medido no sería el de los procesos.

10. **¿Por qué grupos de 16 y no de 64 o de 256?**
    Con grupos grandes, las matrices grupo × fuentes no caben en la caché y varios trabajadores saturan la memoria: con 6 trabajadores, grupos de 64 escalan 1,5 veces y grupos de 16, 3 veces. Con grupos muy chicos (8) pesa el costo de Python de recorrer el árbol muchas veces.

11. **¿Por qué `forkserver` y no `fork`?**
    Con `fork` cada trabajador hereda una copia del proceso padre, incluidos sus arreglos, sin decir nada. Eso escondería cuánto cuesta realmente compartir. Con `forkserver` los trabajadores arrancan limpios, y solo ven lo que se comparte explícitamente.
