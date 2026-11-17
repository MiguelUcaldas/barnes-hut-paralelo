# Fase 5: MPI

## Qué se construyó

Una versión completa en `mpi/` con OpenMPI. Cada proceso tiene su propia memoria: lo único que sabe de los demás es lo que recibe en mensajes.

| Archivo | Contenido |
|---|---|
| `barnes_hut.c` | `main`, paso leapfrog, barrera, energía global, resumen del reparto |
| `particion.c` | cubo raíz global, **particionado de la curva de Morton** (`calcular_cortes`), envío de partículas a su dueño |
| `esencial.c` | cajas de los dominios, **extracción del árbol esencial** para cada proceso, intercambio con `MPI_Alltoallv` |
| `criterio.h` | **criterio de apertura**, el mismo para las fuerzas y para el árbol esencial |
| `arbol.c` / `arbol.h` | árbol local por Morton (como en OpenMP, sin tareas); nodo con hojas vacías, normales y resumidas |
| `fuerzas.c` | recorrido que sirve para el árbol local y para los recibidos; energía local |
| `morton.c` | códigos Morton |
| `entrada_salida.c` | lectura en el proceso 0 y reparto inicial; escritura en el orden del archivo |
| `configuracion.c`, `medicion.c` | opciones; tiempos del proceso más lento, desbalance y CSV |

## Cómo se usa

```bash
make -C mpi
mpirun -np 6 mpi/bin/barnes_hut --entrada datos/grandes/dos_cumulos_200000.txt --pasos 10 --particion costo
mpirun --use-hwthread-cpus -np 12 mpi/bin/barnes_hut ...     # de 7 a 12 procesos
mpirun --oversubscribe -np 16 mpi/bin/barnes_hut ...         # más de 12
```

## Cómo funciona

### El paso completo

```
media patada y deriva (cada proceso, sus partículas)
→ cubo raíz global                   MPI_Allreduce (mínimos y máximos)
→ códigos Morton y cortes de la curva MPI_Allgather (muestras)
→ cada partícula a su dueño           MPI_Alltoall (cantidades) + MPI_Alltoallv (partículas)
→ árbol local
→ cajas de los dominios               MPI_Allgather
→ árbol esencial para cada proceso    MPI_Alltoall (tamaños) + MPI_Alltoallv (nodos y partículas)
→ fuerzas
→ barrera                             MPI_Barrier
→ media patada
```

Todas las comunicaciones son **colectivas**: todos los procesos llaman a la misma función, en el mismo orden. No hay `MPI_Send` / `MPI_Recv` sueltos, así que no puede pasar que dos procesos se queden esperando a que el otro envíe primero (interbloqueo). La condición es que **todos** los procesos ejecuten las mismas colectivas. Por eso funciones como `escribir_aceleraciones_globales` se llaman en todos los procesos, aunque solo el 0 escriba el archivo.

### 1. Particionado por la curva de Morton (`particion.c`)

Cada proceso es dueño de un **tramo contiguo de la curva de Morton**. Un tramo contiguo de la curva es una región compacta del espacio (un conjunto de celdas del octree), así que las partículas de un proceso están juntas y su árbol local es compacto.

**Los cortes** (`calcular_cortes`, la función marcada `PARTICIONADO`) son los códigos donde empieza el tramo de cada proceso. Hay tres formas, con `--particion`:

| Partición | Cómo se cortan los tramos |
|---|---|
| `volumen` | tramos de igual longitud en la curva = igual volumen. No mira dónde están las partículas: es el **reparto ingenuo** |
| `cantidad` | el mismo número de partículas por proceso |
| `costo` (por defecto) | el mismo número de **interacciones del paso anterior** por proceso |

Para `cantidad` y `costo`, ningún proceso conoce todas las partículas. Cada uno:

1. ordena las suyas por código;
2. toma 256 **muestras** repartidas a igual peso (`tomar_muestras`). Cada muestra es un código y el peso que representa: 1 por partícula para `cantidad`, sus interacciones para `costo`;
3. con `MPI_Allgather`, todos reciben todas las muestras (256 × P, unos pocos KB);
4. todos ordenan las mismas muestras y cortan donde el peso acumulado llega a 1/P, 2/P, …, así que **todos calculan exactamente los mismos cortes** sin que nadie tenga que decidir por los demás.

**El envío** (`intercambiar_particulas`). Como las partículas locales están ordenadas, las que van a cada destino están seguidas. Primero un `MPI_Alltoall` le avisa a cada proceso cuántos bytes va a recibir de cada uno, porque sin eso no podría reservar memoria. Después un solo `MPI_Alltoallv` mueve todas las partículas. La estructura `Particula` lleva todo el estado (posición, velocidad, masa, costo, identificador) para que una partícula viaje en un solo bloque.

### 2. Árbol local

Cada proceso construye el octree de **sus** partículas con la construcción por Morton de la fase 3, sin tareas. Lo importante es que **el cubo raíz es el global**, igual en todos los procesos: una celda del árbol del proceso A y una del proceso B con el mismo código son exactamente el mismo cubo del espacio.

### 3. Árbol esencial local (`esencial.c`)

Para calcular la fuerza sobre sus partículas, un proceso necesita información de las partículas de los demás. Pero el enunciado pide **no replicar el árbol completo**. La idea:

- Las partículas lejanas se pueden resumir: el proceso solo necesita la masa y el centro de masa de una celda lejana, no sus partículas.
- Lo que está "lejos" se decide con el criterio θ.

Cada proceso publica la **caja de su dominio** (mínimo y máximo de sus partículas, con `MPI_Allgather`). Después, para cada otro proceso r, recorre **su propio** árbol con el criterio evaluado contra **la caja entera de r** (`extraer_arbol_podado`):

| El nodo… | Se envía… |
|---|---|
| está lejos de toda la caja de r | **resumido**: masa y centro de masa, sin hijos (hoja "resumida") |
| es una hoja cercana | con sus partículas (posición y masa) |
| es interno y está cerca | el nodo, y se siguen revisando sus 8 hijos |

El árbol podado se arma en un arreglo que también sirve de cola: los nodos se agregan al final y se procesan en orden. Así los 8 hijos de cada nodo quedan seguidos, como en el original, y el árbol recibido se recorre con el **mismo código** que el árbol local.

**Por qué alcanza con lo que se envía.** El criterio (`criterio.h`) mide la distancia desde el centro de masa de la celda hasta el punto más cercano de una caja:

- al calcular fuerzas, la caja es la partícula (mínimo = máximo);
- al extraer el árbol esencial, la caja es todo el dominio del otro proceso.

Una partícula de r está dentro de la caja de r, así que nunca está más cerca de una celda que la caja. Si una celda está lejos de la caja, también está lejos de cada partícula de r. Por eso r nunca necesita abrir una celda que le llegó resumida. Como es **la misma función** para los dos usos, si en la demostración se cambia el criterio, la garantía se mantiene sola.

**El intercambio** (`enviar_a_todos`): se junta en un solo buffer lo que va a cada destino, un `MPI_Alltoall` avisa los tamaños, y un `MPI_Alltoallv` envía. Se hace dos veces: una para los nodos y otra para las partículas de las hojas. Los índices dentro de cada árbol recibido (`primer_hijo`, `primera_particula`) son relativos a su propio inicio, así que el árbol se usa directamente sobre el buffer, sin traducir nada.

### 4. Fuerzas (`fuerzas.c`)

Cada partícula recorre su árbol local (saltándose a sí misma) y luego el árbol esencial recibido de cada uno de los otros procesos. El recorrido es uno solo:

1. Si la celda está vacía, se ignora.
2. Si está lejos según el criterio, o es una hoja resumida, actúa como una masa puntual.
3. Si es una hoja normal, se interactúa con cada partícula.
4. Si no, se baja a los 8 hijos.

### 5. La barrera y el desbalance

Después de las fuerzas hay un `MPI_Barrier`. **En MPI no es necesario para que el resultado sea correcto**: cada proceso solo modifica sus propias partículas, y lo que recibió de los demás es una copia. Se pone para medir bien. Sin la barrera, un proceso rápido llegaría antes al `MPI_Allreduce` del paso siguiente y esperaría ahí; esa espera quedaría contada como "orden y reparto". Con la barrera, la espera queda en la fase de fuerzas, que es donde se produjo.

El desbalance se mide con el tiempo de cálculo propio de cada proceso (antes de la barrera): máximo / promedio, reunidos con `MPI_Reduce`.

### 6. Entrada y salida

El proceso 0 lee el archivo y reparte bloques iguales de filas (`MPI_Scatterv`). En el primer cálculo de fuerzas se reparten de nuevo por la curva de Morton. Para escribir, todos envían sus partículas al 0 (`MPI_Gatherv`), que las ordena por identificador.

Si el archivo no existe, solo el proceso 0 se entera, y los demás quedarían esperando en el `MPI_Bcast`: por eso se usa `MPI_Abort`, que termina a todos.

## Un problema de rendimiento que se encontró: `fmax`

La primera versión del criterio usaba `fmax` para la distancia a la caja. Con 1 proceso, el cálculo de fuerzas tardaba 3,78 s por paso, contra 1,87 s de la versión OpenMP con 1 hilo, con el mismo árbol.

`fmax` debe tratar los NaN de una forma especial, y por eso gcc no la convierte en una sola instrucción. Se reemplazó por `a > b ? a : b` (`distancia_a_intervalo`), que corresponde a la instrucción `maxsd` del procesador. Así bajó a 2,22 s, a 5 % de una versión que solo sirve para puntos.

## Resultados (WSL, Ryzen 5 5500U: 6 núcleos, 12 hilos)

### Correctitud

| | θ = 0 | θ = 0,5, 1 proceso | θ = 0,5, 16 procesos |
|---|---|---|---|
| uniforme | 1,4e-10 | 0,48 % | 0,49 % |
| dos_cumulos | 7,3e-11 | 0,41 % | 0,44 % |
| concentrada | 9,1e-11 | 0,36 % | 0,41 % |

- **Con θ = 0 el resultado es la suma directa con 1, 3, 8 y 16 procesos.** Con 1 proceso, el error con θ = 0,5 es idéntico al de C.
- **Con más procesos, el error con θ = 0,5 sube un poco** (de 0,48 % a 0,49 % en uniforme). Una celda que en el árbol global tendría partículas de dos procesos aquí son dos celdas, una por proceso. Cada parte tiene el mismo tamaño pero su masa está cargada hacia un lado, y aproximarla por su centro de masa es un poco peor.
- **Energía** en 100 pasos con 6 procesos, según `verificar.py`: 2,9e-6 (uniforme), 3,7e-4 (dos_cumulos), 1,2e-4 (concentrada). Todas pasan.
- **Casos límite**: con 4 partículas y 6 procesos, varios procesos quedan sin partículas y todo sigue funcionando. Las partículas repetidas y todas en un punto también funcionan.
- **`make sanitizado`** (UBSan) con 4 y 16 procesos y las tres particiones: sin avisos.

**AddressSanitizer no se pudo usar.** En este sistema, cualquier programa MPI compilado con `-fsanitize=address` falla dentro de `MPI_Init`, incluso uno vacío, así que es una incompatibilidad entre AddressSanitizer y OpenMPI 5. Para revisar la memoria en MPI haría falta valgrind (`sudo apt install valgrind`).

### Escalado (N = 200.000, partición por costo, ms por paso, una corrida de 3 pasos)

| Procesos | uniforme | aceleración | árbol esencial | dos_cumulos | aceleración | árbol esencial |
|---|---|---|---|---|---|---|
| 1 | 2336 | 1,00 | 3 | 5103 | 1,00 | 3 |
| 2 | 1360 | 1,72 | 13 | 3077 | 1,66 | 48 |
| 4 | 925 | 2,52 | 38 | 1981 | 2,58 | 113 |
| 6 | 811 | 2,88 | 94 | 1615 | 3,16 | 139 |
| 8 | 686 | 3,40 | 83 | 1457 | 3,50 | 165 |
| 12 | 686 | 3,41 | 128 | 1415 | 3,61 | 211 |
| 16 | 698 | 3,35 | 105 | 1565 | 3,26 | 229 |

**Cuánto árbol recibe cada proceso** (uniforme): el árbol completo tiene unos 751.000 nodos.

| Procesos | árbol local | árbol esencial recibido (promedio / máximo) |
|---|---|---|
| 2 | 375.000 | 58.000 / 81.000 |
| 6 | 125.000 | 200.000 / 332.000 |
| 12 | 63.000 | 145.000 / 394.000 |

Nadie recibe el árbol completo: con 12 procesos, un proceso tiene en promedio 63.000 + 145.000 = 208.000 nodos (28 % del total).

**Por qué escala menos que OpenMP** (con 12, OpenMP tarda 511 ms en uniforme y 960 ms en dos_cumulos; MPI, 686 y 1415):

- **El árbol esencial cuesta, y crece con P.** Con 12 procesos son 128–211 ms por paso (19 % y 15 %). Cada proceso recorre su árbol una vez por cada otro proceso y envía el resultado. Con P procesos hay P·(P−1) árboles podados por paso.
- **Más trabajo en las fuerzas.** Cada partícula recorre P árboles (el local y P−1 recibidos) en lugar de uno. Las celdas de la frontera entre dominios llegan partidas en varias, una por proceso, y esto agrega interacciones.
- **Los procesos no comparten memoria.** El árbol esencial se copia en cada proceso. Con 12 procesos hay 12 × 208.000 nodos en memoria, contra 751.000 con OpenMP.
- **De 8 a 12**, igual que en las otras versiones, los hilos SMT aportan poco. **Con 16** hay más procesos que hilos de hardware: el desbalance sube a 1,20 y el tiempo empeora.

### Comparación de particiones (N = 200.000, ms por paso)

| Distribución | Procesos | volumen | cantidad | costo |
|---|---|---|---|---|
| dos_cumulos | 6 | 2411 (desb. 1,69) | 1654 (1,08) | 1658 (1,06) |
| dos_cumulos | 12 | 2216 (1,77) | 1398 (1,07) | 1399 (1,06) |
| concentrada | 6 | 6074 (4,31) | 1751 (1,17) | 1691 (1,08) |
| concentrada | 12 | 9271 (**8,66**) | 1601 (1,21) | 1444 (1,10) |

Partículas por proceso con `concentrada` y 12 procesos:

- `volumen`: de 524 a 125.945. Un proceso se queda con el 63 % de las partículas y los demás lo esperan. **Con 12 procesos es más lento que con 6**: agregar procesos solo agrega procesos vacíos.
- `cantidad`: de 16.564 a 16.734, pero las del cúmulo cuestan más, así que el desbalance queda en 1,21.
- `costo`: de 12.785 a 26.855 partículas. Reparte **trabajo** y no partículas: los procesos del cúmulo reciben menos. Es el más rápido.

## Qué se ve en htop

P procesos `barnes_hut`, cada uno en un núcleo, en 100 % durante las fuerzas. Con `--particion volumen` y `concentrada`, uno o dos procesos siguen en 100 % mientras los demás bajan. El que espera en `MPI_Barrier` también aparece ocupado, porque OpenMPI espera revisando activamente si llegó el mensaje, sin dormir. Por eso la columna de CPU de htop no distingue un proceso que trabaja de uno que espera: hay que mirar el desbalance que imprime el programa. Con `--use-hwthread-cpus`, los procesos se reparten también en los hilos SMT.

## Dónde se hacen las modificaciones en vivo

| Pedido | Dónde | Qué cambiar |
|---|---|---|
| Cambiar la estrategia de particionado | `particion.c`, `calcular_cortes` | por ejemplo, otro peso en `peso`, o cortes fijos |
| Particionar por bisección recursiva en vez de Morton | `particion.c` | otra función de cortes que divida por la mediana en el eje más largo; el resto (envío, árbol esencial) no cambia, porque solo usa las cajas de los dominios |
| Criterio θ | `criterio.h`, `celda_suficientemente_lejana` | se aplica solo a las fuerzas y al árbol esencial |
| Muestras para los cortes | `particion.c`, `MUESTRAS_POR_PROCESO` | más muestras: cortes más exactos y un `Allgather` más grande |

Ojo con la bisección: el envío actual supone que cada proceso recibe un tramo de la curva (`cortes` crecientes). Con bisección, el destino de cada partícula se calcularía con su posición y no con su código.

## Decisiones de diseño

- **Curva de Morton y no bisección recursiva.** El orden de Morton ya está hecho (fases 3 y 4), un tramo de la curva es una región compacta, y cortar la curva se reduce a elegir P − 1 números. La bisección da dominios en forma de caja, que son más compactos (cajas más chicas, árboles esenciales más pequeños), pero necesita varias rondas de comunicación por nivel.
- **Cortes por muestreo.** Nadie necesita conocer todas las partículas: alcanzan 256 muestras por proceso.
- **Colectivas en lugar de mensajes punto a punto.** Son más simples de razonar y no pueden interbloquearse si todos las llaman en el mismo orden.
- **Un solo criterio para las fuerzas y para el árbol esencial**, con la distancia a una caja. Así es imposible que el árbol esencial quede incompleto después de cambiar el criterio.
- **Partición por costo como valor por defecto**: fue la más rápida en las distribuciones con cúmulos.

## Preguntas tipo sustentación

1. **¿Qué es el árbol esencial y por qué no se envía el árbol completo?**
   Es la parte del árbol de un proceso que otro proceso necesita: las celdas lejanas resumidas en su masa y centro de masa, y las cercanas con detalle. Enviar el árbol completo a cada proceso costaría P veces su tamaño en comunicación y memoria. Con 12 procesos, cada uno recibe en promedio el 19 % del árbol.

2. **¿Cómo sabe un proceso qué parte de su árbol necesita otro?**
   Por la caja del dominio del otro (`MPI_Allgather`). Recorre su árbol con el criterio evaluado contra la caja entera: lo que está lejos de toda la caja se envía resumido; lo demás, con detalle.

3. **¿Por qué es seguro enviar una celda resumida?**
   Porque toda partícula del otro proceso está dentro de su caja, y entonces su distancia a la celda es mayor o igual que la distancia de la caja. Si la celda cumple el criterio para la caja, lo cumple para cada partícula, y nadie va a necesitar abrirla.

4. **¿Cómo se decide qué partículas tiene cada proceso?**
   Cada uno ordena sus partículas por código Morton y aporta 256 muestras con su peso (`MPI_Allgather`). Todos ordenan las mismas muestras y cortan donde el peso acumulado llega a 1/P, 2/P…, así que todos obtienen los mismos cortes. Luego cada partícula se envía a su dueño con `MPI_Alltoallv`.

5. **¿Por qué la partición por volumen es tan mala con la distribución concentrada?**
   Corta el espacio en partes iguales sin mirar las partículas. Con el 90 % de las partículas en un cúmulo, un proceso se queda con 126.000 de 200.000 y el resto espera. Con 12 procesos el desbalance es 8,66 y es más lento que con 6.

6. **¿Puede haber interbloqueo en este programa? ¿Cuándo?**
   Todas las comunicaciones son colectivas. Habría interbloqueo si un proceso no llama a una colectiva que los demás sí llaman. Por ejemplo, escribir `if (rango == 0) escribir_aceleraciones_globales(...)`: el 0 esperaría en el `MPI_Gatherv` a procesos que nunca llegan. Por eso esa función la llaman todos.

7. **¿Por qué hay un `MPI_Alltoall` antes de cada `MPI_Alltoallv`?**
   Porque el receptor tiene que saber cuántos bytes va a recibir de cada proceso para reservar memoria y calcular dónde poner cada bloque. El primero intercambia los tamaños y el segundo los datos.

8. **¿Es necesaria la barrera después de las fuerzas?**
   Para la correctitud, no: cada proceso solo modifica sus partículas y los datos remotos son copias. Se pone para medir: sin ella, la espera de los procesos rápidos se contaría en la fase siguiente.

9. **¿Por qué MPI escala menos que OpenMP en el mismo equipo?**
   Por el costo del árbol esencial (extraerlo y enviarlo crece con el número de procesos), porque cada partícula recorre P árboles y las celdas de frontera llegan partidas, y porque cada proceso tiene su propia copia de lo que recibe. OpenMP comparte un solo árbol en memoria.

10. **¿Qué pasa si el archivo de entrada no existe?**
    Solo el proceso 0 lo lee y se da cuenta. Si terminara solo él, los demás quedarían esperando para siempre en el `MPI_Bcast`. Por eso se usa `MPI_Abort`, que termina todos los procesos.

11. **¿Por qué con más procesos el error de fuerzas sube un poco?**
    Una celda del espacio que tiene partículas de varios procesos llega como varias celdas, una por proceso. Cada parte tiene el tamaño de la celda completa, pero su masa está cargada hacia un lado. Aproximarla por su centro de masa es un poco menos preciso que aproximar la celda entera.

12. **¿Por qué con 16 procesos va más lento que con 12?**
    Hay 12 hilos de hardware. Con 16 procesos, el sistema operativo los turna. Un proceso que pierde su turno retrasa a todos, porque los demás lo esperan en la siguiente colectiva. El desbalance sube a 1,20.
