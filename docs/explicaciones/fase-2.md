# Fase 2: Barnes-Hut secuencial en C

## Qué se construyó

Un programa completo en `c_openmp/` que lee las partículas, arma el octree en cada paso, calcula las fuerzas con el criterio θ, integra con leapfrog y mide el tiempo de cada fase. Es la base de la versión OpenMP (fase 3, en esta misma carpeta) y de la versión MPI (fase 5, que parte de una copia del árbol).

| Archivo | Contenido |
|---|---|
| `barnes_hut.c` | `main` y el paso leapfrog completo (`avanzar_un_paso`, `media_patada`, `deriva`) |
| `arbol.c` / `arbol.h` | `Nodo`, `Arbol`, cubo raíz, inserción, subdivisión y masas |
| `fuerzas.c` / `fuerzas.h` | **criterio de apertura**, interacción con suavizado, recorrido del árbol, energías |
| `particulas.c` / `particulas.h` | `Vector3`, `Particulas`, lectura y escritura de archivos |
| `configuracion.c` / `configuracion.h` | opciones de línea de comandos |
| `medicion.c` / `medicion.h` | cronómetro, impresión de tiempos y fila del CSV |
| `memoria.h` | `calloc` / `realloc` que terminan con un mensaje si falta memoria |
| `Makefile` | `make` (normal, `-O2`) y `make sanitizado` (detector de errores de memoria) |

## Cómo se usa

```bash
cd c_openmp && make && cd ..
c_openmp/bin/barnes_hut --entrada datos/uniforme_2000.txt --pasos 100 --theta 0.5 --epsilon 0.01 --dt 0.001
```

Opciones: `--entrada`, `--pasos` (100), `--dt` (0.001), `--theta` (0.5), `--epsilon` (0.01), `--cada-energia K`, `--guardar-aceleraciones f`, `--guardar-final f`, `--distribucion nombre`, `--csv f` (`resultados/mediciones.csv`) y `--ayuda`. El programa se corre desde la raíz del repositorio para que la ruta del CSV funcione.

Verificación:

```bash
c_openmp/bin/barnes_hut --entrada datos/uniforme_2000.txt --pasos 100 \
    --guardar-aceleraciones /tmp/acel.txt --guardar-final /tmp/final.txt
python comun/verificar.py fuerzas --particulas datos/uniforme_2000.txt --aceleraciones /tmp/acel.txt --epsilon 0.01
python comun/verificar.py energia --inicial datos/uniforme_2000.txt --final /tmp/final.txt --epsilon 0.01
```

## Cómo funciona

### El árbol en un arreglo plano

Cada `Nodo` guarda el centro geométrico del cubo, la mitad del lado, la masa, el centro de masa, `primer_hijo` y `particula`. Los nodos están en un solo arreglo (`arbol->nodos`) y se referencian por índice. El nodo 0 es la raíz.

Cuando un nodo se divide, sus **8 hijos se crean seguidos**: el hijo k está en `primer_hijo + k`. Así basta un entero para llegar a los ocho.

El número k de cada hijo tiene un bit por eje: bit 0 = x, bit 1 = y, bit 2 = z (1 = lado positivo). Por ejemplo, el hijo 5 = 101 en binario está del lado positivo en x y en z, y del lado negativo en y. La función `octante` calcula ese número para una posición con tres comparaciones.

¿Por qué un arreglo y no un `malloc` por nodo?

- Se reserva una vez y se reutiliza en todos los pasos. Reconstruir es solo poner `num_nodos = 0`.
- Los nodos quedan juntos en memoria, lo que aprovecha mejor la caché.
- Un arreglo se puede enviar por MPI o copiar a la GPU tal cual. Una estructura con punteros, no.

### Construcción (`construir_arbol`)

1. **Cubo raíz** (`calcular_cubo_raiz`): mínimo y máximo de las posiciones en cada eje. El centro es el punto medio y el lado es la mayor de las tres extensiones. Se recalcula en cada paso.
2. **Inserción** (`insertar`): cada partícula baja desde la raíz eligiendo en cada nodo interno el hijo de su octante.
   - Si llega a una hoja vacía, se queda ahí.
   - Si la hoja está ocupada, la hoja se divide en 8. La partícula que estaba baja a su octante y la nueva sigue bajando. Si las dos caen en el mismo octante, se vuelve a dividir, hasta que quedan separadas.
3. **Masas** (`calcular_masas`): en una hoja, la masa y el centro de masa son los de su partícula. En un nodo interno, la masa es la suma de las de sus hijos y el centro de masa es el promedio de los centros de masa de los hijos, ponderado por su masa.

**El truco del recorrido hacia atrás.** Un hijo siempre se crea después que su padre, así que tiene un índice mayor. Si se recorre el arreglo del último nodo al primero, cuando se llega a un padre todos sus hijos ya tienen su masa calculada. Es un recorrido de abajo hacia arriba sin recursión.

**Cuidado con `realloc`.** Si el arreglo se llena, `nuevo_nodo` duplica su capacidad con `realloc`, y eso puede mover todo el arreglo a otra dirección. Cualquier puntero `Nodo *` guardado de antes queda apuntando a memoria liberada. Por eso `subdividir` copia el centro y el lado del padre antes de crear hijos, y `insertar` vuelve a tomar el puntero (`&arbol->nodos[actual]`) después de subdividir. Es un error clásico que el modo sanitizado detecta.

**Partículas en el mismo punto.** Si dos partículas tienen exactamente la misma posición, dividir nunca las separa y el árbol crecería sin fin. Al llegar a `PROFUNDIDAD_MAXIMA` (60), la hoja guarda varias partículas encadenadas con el arreglo `siguiente`. A profundidad 60, el cubo mide 2⁻⁶⁰ ≈ 1e-18 veces el cubo raíz, así que solo llegan ahí partículas prácticamente idénticas.

### Fuerzas (`calcular_aceleraciones`)

Para cada partícula i se recorre el árbol desde la raíz (`recorrer_arbol`, recursivo):

- **Hoja**: interacción exacta con su partícula, saltándose a sí misma (`j != i`).
- **Nodo interno**: se calcula la distancia d de i al centro de masa del nodo.
  - Si `tamaño / d < θ` (`celda_suficientemente_lejana`), toda la celda actúa como una partícula de masa `nodo->masa` en su centro de masa.
  - Si no, se baja a los 8 hijos.

La interacción (`sumar_atraccion`) usa el suavizado: `a += m·d / (|d|² + ε²)^(3/2)` y `φ -= m / √(|d|² + ε²)`. Calcula una sola raíz y de ella saca 1/r y 1/r³.

**Criterio con cuadrados.** `tamaño/d < θ` equivale a `tamaño² < θ²·d²`, porque todo es positivo. Así no hace falta una raíz por cada celda visitada, solo por cada interacción.

**Por qué el ciclo será fácil de paralelizar.** Cada partícula solo *lee* el árbol y las posiciones, y solo *escribe* su propia aceleración y su potencial. No hay dos iteraciones que escriban lo mismo, así que no hay condición de carrera. Esto es lo que se aprovecha en la fase 3.

### Energía

Al recorrer el árbol, además de la aceleración se acumula el potencial φ_i de cada partícula. La energía potencial es `U = ½ Σ m_i φ_i`. El ½ está porque cada par aparece en φ_i y en φ_j.

Así la energía no cuesta otro recorrido. Sale aproximada con el mismo θ, y con 200.000 partículas una suma directa tardaría casi una hora. Para la verificación exacta con N pequeño está `verificar.py`.

### El paso de tiempo (`avanzar_un_paso`)

```
media_patada → deriva → construir_arbol → calcular_aceleraciones → media_patada
```

Antes del ciclo se calculan una vez las aceleraciones del estado inicial: la primera media patada las necesita. Ese cálculo inicial se informa aparte y no entra en los tiempos por paso.

**Foto congelada.** Todas las aceleraciones de un paso se calculan con el mismo árbol, armado con las mismas posiciones. Nadie se mueve hasta que se terminó de calcular la fuerza sobre todos. En la versión secuencial esto se cumple solo; en las paralelas exigirá una barrera entre fases.

### Medición

`segundos_actuales` usa `clock_gettime(CLOCK_MONOTONIC)`, que no salta si el reloj del sistema se ajusta durante la corrida.

Cada paso mide construcción, fuerzas, integración (las dos medias patadas y la deriva) y total. Ordenamiento queda en 0 porque esta versión no ordena. Al final se imprimen los acumulados y una fila en el CSV. La columna `version` es `c_secuencial`.

## Dónde se hacen las modificaciones en vivo

| Pedido | Dónde | Qué cambiar |
|---|---|---|
| Cambiar el criterio θ | `fuerzas.c`, `celda_suficientemente_lejana` | Por ejemplo, medir la distancia al centro geométrico en vez de al centro de masa: calcular `distancia2` con `nodo->centro` en `recorrer_arbol`. O usar la diagonal del cubo (`tamano = 2·√3·mitad_lado`) |
| Cambiar θ, ε, Δt o pasos | línea de comandos | `--theta`, `--epsilon`, `--dt`, `--pasos` |
| Más profundidad o celdas con varias partículas | `arbol.h`, `PROFUNDIDAD_MAXIMA`; `arbol.c`, `insertar` | — |
| Cambiar el integrador | `barnes_hut.c`, `avanzar_un_paso` | — |

## Resultados (WSL, Ryzen 5 5500U, un hilo)

### Error de fuerzas contra la suma directa (N = 2.000, ε = 0,01)

| θ | uniforme (prom. / máx.) | dos_cumulos | concentrada |
|---|---|---|---|
| 0 | 1,4e-10 / 4,2e-10 | 7,3e-11 / 4,8e-10 | 9,1e-11 / 4,6e-10 |
| 0,3 | 0,12 % / 2,5 % | 0,09 % / 0,9 % | 0,08 % / 0,8 % |
| 0,5 | 0,48 % / 5,2 % | 0,41 % / 2,8 % | 0,36 % / 5,2 % |
| 0,7 | 1,24 % / 20 % | 0,98 % / 7,9 % | 0,77 % / 4,2 % |
| 1,0 | 3,4 % / 37 % | 2,4 % / 28 % | 1,7 % / 11 % |

Con θ = 0 el error es solo el redondeo del archivo: el recorrido baja siempre hasta las hojas y equivale a la suma directa. El error crece con θ, como se espera.

### Energía (N = 2.000, θ = 0,5, Δt = 0,001, 100 pasos)

| Distribución | Error según el programa (árbol) | Error según `verificar.py` (exacto) | Referencia con fuerzas exactas (fase 1) |
|---|---|---|---|
| uniforme | 2,1e-5 | 6,3e-6 | 3,4e-8 |
| dos_cumulos | 1,4e-4 | 2,4e-4 | 4,8e-5 |
| concentrada | 1,4e-4 | 1,2e-4 | 5,9e-5 |

Las tres pasan la tolerancia de 1e-3. El error es mayor que con fuerzas exactas porque la aproximación del árbol no es exactamente conservativa, pero sigue siendo pequeño.

### Tiempo por paso (θ = 0,5, promedio de 5 pasos)

| Distribución | N | construcción | fuerzas | total | nodos | profundidad |
|---|---|---|---|---|---|---|
| uniforme | 20.000 | 4,4 ms | 295 ms | 0,30 s | 75.161 | 12 |
| uniforme | 200.000 | 81 ms | 6,58 s | 6,67 s | 750.905 | 12 |
| dos_cumulos | 20.000 | 4,5 ms | 611 ms | 0,62 s | 78.233 | 13 |
| dos_cumulos | 200.000 | 110 ms | 11,3 s | 11,4 s | 771.705 | 15 |
| concentrada | 20.000 | 5,7 ms | 623 ms | 0,63 s | 77.497 | 13 |
| concentrada | 200.000 | 107 ms | 10,4 s | 10,5 s | 774.849 | 15 |

Observaciones:

- **El 99 % del tiempo está en las fuerzas.** Ahí está casi toda la ganancia posible al paralelizar. La construcción es ~1 %, pero en la versión paralela pesará más: si las fuerzas se dividen entre 8 hilos y la construcción no, la construcción pasa a ser cerca del 8 % (ley de Amdahl). Por eso el enunciado pide construirla en paralelo.
- **Con 10 veces más partículas, el tiempo sube unas 22 veces** (uniforme): N log N predice ~12. El resto se explica por la caché. Las partículas están en orden aleatorio, así que dos partículas seguidas recorren partes distintas del árbol y casi todo acceso a un nodo va a la RAM. Con 750.000 nodos de 72 bytes (54 MB), el árbol no cabe en la caché L3 del procesador (4 MB según `lscpu` en WSL). Ordenar las partículas por Morton (fase 3) hace que partículas consecutivas recorran casi los mismos nodos.
- **Una partícula de un cúmulo cuesta el doble que una aislada.** En las zonas densas hay más celdas cerca, y las celdas cercanas no cumplen el criterio y se abren. Por eso repartir las partículas en partes iguales entre hilos no reparte el trabajo en partes iguales.

### Errores de memoria

Con `make sanitizado` (AddressSanitizer y UndefinedBehaviorSanitizer) se corrieron 20 pasos con dos cúmulos y tres casos límite: partículas repetidas, dos partículas a 1e-12 (profundidad 40, el arreglo de nodos crece con `realloc` de 20 a 321 nodos) y todas en un punto. No hubo ningún aviso, y con θ = 0 los tres coinciden con la suma directa.

## Decisiones de diseño

- **Inserción de arriba hacia abajo**, como la describe el enunciado. Es la forma más fácil de explicar. En la fase 3 la construcción paralela usará candados o Morton.
- **Hijos contiguos** (`primer_hijo + k`) en vez de un arreglo `hijos[8]`. Ocupa menos memoria por nodo (4 bytes en vez de 32) y el árbol queda listo para aplanarse en CUDA. El costo es crear 8 hijos aunque algunos queden vacíos.
- **Recorrido recursivo** para las fuerzas. Es el más claro de leer y la profundidad (≤ 60) no es un problema para la pila. En CUDA se cambiará por una pila explícita.
- **Distancia al centro de masa** en el criterio, como dice el enunciado. Con θ ≤ 1/√3 ≈ 0,577, una celda que contiene a la partícula nunca se aproxima. La partícula y el centro de masa están dentro del mismo cubo, así que d ≤ √3·tamaño, y entonces tamaño/d ≥ 1/√3 > θ. Con θ más grande sí puede pasar, y es una de las causas de que el error máximo crezca rápido con θ = 0,7 y 1.
- **Energía con el mismo árbol**, sin otro recorrido. Es aproximada, pero sirve para N grande. La verificación exacta está en `verificar.py`.
- **Se compila con `-O2`.** `-O3 -march=native` puede ser más rápido, pero el binario dependería del procesador. Se puede probar en la fase 7.

## Preguntas tipo sustentación

1. **¿Qué pasa con el resultado si θ = 0? ¿Y con el tiempo?**
   Ninguna celda cumple `tamaño² < 0`, así que todo se abre hasta las hojas. El resultado es la suma directa: error de 1e-10 contra la referencia. El tiempo es O(N²) y además más lento que una suma directa simple, porque se recorre el árbol completo para cada partícula.

2. **¿Por qué el criterio se escribe con cuadrados?**
   `s/d < θ` equivale a `s² < θ²d²` porque todo es positivo. Así se evita calcular una raíz cuadrada por cada celda visitada.

3. **¿Cómo se calcula la masa de un nodo interno y por qué se recorre el arreglo al revés?**
   Es la suma de las masas de sus 8 hijos, y su centro de masa es el promedio de los centros de masa de los hijos ponderado por su masa. Los hijos se crean después que el padre y tienen índices mayores. Recorriendo de atrás hacia adelante, los hijos ya están listos cuando se llega al padre, sin recursión.

4. **¿Qué error se cometería guardando un `Nodo *` antes de llamar a `subdividir`?**
   `subdividir` crea nodos y puede llamar a `realloc`, que mueve el arreglo. El puntero viejo apuntaría a memoria liberada: el programa puede funcionar con N pequeño y fallar o dar valores basura con N grande. Por eso se usan índices.

5. **¿Qué pasa si dos partículas están en la misma posición?**
   Dividir nunca las separaría. Al llegar a profundidad 60 la hoja guarda las dos, encadenadas con `siguiente`. La fuerza entre ellas la limita el suavizado ε: su distancia es 0, pero el denominador es ε³ y no 0.

6. **¿Por qué el cálculo de fuerzas se puede paralelizar sin candados?**
   Cada iteración del ciclo de partículas lee datos compartidos (árbol y posiciones) que nadie modifica en esa fase, y escribe solo `aceleraciones[i]` y `potenciales[i]`, que son distintos para cada i. No hay escrituras compartidas.

7. **Si la construcción es el 1 % del tiempo, ¿por qué el enunciado exige construir el árbol en paralelo?**
   Por la ley de Amdahl: si las fuerzas se aceleran 8 veces y la construcción no, la parte secuencial pasa del 1 % a cerca del 8 % del total. Con más hilos sigue creciendo y limita la aceleración máxima.

8. **¿Por qué una distribución con cúmulos tarda el doble que la uniforme con el mismo N?**
   Las partículas de un cúmulo tienen muchas celdas cerca, y las celdas cercanas no cumplen el criterio, así que se abren. Cada partícula visita más nodos. Además el árbol es más profundo (15 contra 12 niveles).

9. **¿Por qué al pasar de 20.000 a 200.000 el tiempo sube 22 veces y no unas 12, como predice N log N?**
   Por la memoria. El árbol de 200.000 partículas ocupa 54 MB y no cabe en la caché (4 MB de L3). Las partículas están en orden aleatorio, así que cada una recorre partes distintas del árbol y casi todo acceso va a la RAM. Ordenarlas por posición (Morton) lo mejora.

10. **¿Qué significa "foto congelada" y cómo se cumple aquí?**
    Que todas las aceleraciones de un paso se calculan con las mismas posiciones. Aquí el árbol se construye, luego se calculan todas las aceleraciones y solo después se mueven las partículas. En las versiones paralelas habrá que poner una barrera para que ningún hilo mueva partículas mientras otro aún calcula fuerzas.
