# Fase 1: herramientas comunes

## Qué se construyó

En `comun/` quedan las herramientas que usan las cuatro versiones. Ninguna hace parte de la simulación paralela.

| Archivo | Para qué sirve |
|---|---|
| `formato.py` | Lee y escribe los archivos de partículas (7 columnas) y de aceleraciones (3 columnas) |
| `generar_particulas.py` | Crea las distribuciones `uniforme`, `dos_cumulos` y `concentrada` con semilla fija |
| `referencia_directa.py` | Suma directa O(N²): aceleraciones exactas, energía K + U y una simulación leapfrog con fuerzas exactas |
| `verificar.py` | Compara la salida de cualquier versión con la suma directa: error de fuerzas y error de energía |

Además, `scripts/generar_datos.sh` genera todos los archivos: los de 2.000 partículas en `datos/`, que se suben a git, y los de 20.000 a 200.000 en `datos/grandes/`, que no se suben.

## Cómo se usa

```bash
# Generar
python comun/generar_particulas.py --distribucion dos_cumulos --n 2000
bash scripts/generar_datos.sh .venv/bin/python

# Referencia exacta
python comun/referencia_directa.py energia --particulas datos/uniforme_2000.txt --epsilon 0.01
python comun/referencia_directa.py simular --particulas datos/uniforme_2000.txt --epsilon 0.01 --dt 0.001 --pasos 100

# Verificar una versión (cuando exista)
python comun/verificar.py fuerzas --particulas datos/uniforme_2000.txt --aceleraciones acel.txt --epsilon 0.01
python comun/verificar.py energia --inicial datos/uniforme_2000.txt --final final.txt --epsilon 0.01
```

`verificar.py` termina con código de salida 0 si la versión pasa y 1 si no. Así los scripts de bash pueden detenerse si una versión falla.

## Lo que deben escribir las versiones para poder verificarlas

Cada versión (C, OpenMP, MPI, Python y CUDA) va a aceptar dos opciones:

- `--guardar-aceleraciones archivo`: escribe las aceleraciones del estado **inicial** (las que se calculan antes del primer paso), una fila `ax ay az` por partícula, en el orden del archivo de entrada.
- `--guardar-final archivo`: escribe el estado final en el mismo formato de 7 columnas de la entrada.

Con eso, `verificar.py` revisa cualquier versión de la misma forma, sin importar el lenguaje.

## Cómo funciona

### Las distribuciones

- **uniforme**: posiciones al azar en el cubo unitario, velocidades cero y masas 1/N. Es la que pide el enunciado para desarrollar.
- **dos_cumulos**: dos esferas de Plummer de 1.000 o 10.000 partículas cada una. Están separadas 0,4 en x y se acercan con rapidez relativa 0,5.
- **concentrada**: el 90 % de las partículas en una esfera de Plummer de radio de escala 0,05 centrada en (0,3; 0,3; 0,3). El otro 10 % queda repartido al azar por el cubo.

Una **esfera de Plummer** es un modelo clásico de cúmulo estelar. La densidad es alta en el centro y baja hacia afuera. Las posiciones se sacan invirtiendo la fórmula de la masa acumulada. Las velocidades se sacan con el método de Aarseth, Hénon y Wielen (1974): cada partícula recibe una fracción de la velocidad de escape a su radio. Esa fracción sigue la curva q²(1 − q²)^(7/2) y se muestrea por rechazo. Con estas velocidades el cúmulo está en equilibrio: no colapsa ni se dispersa.

Todas las filas se barajan al final. Así ninguna versión puede depender, sin darse cuenta, de que partículas vecinas en el archivo estén cerca en el espacio. Ordenar las partículas por posición (Morton) es trabajo de cada versión.

### Por qué estas distribuciones castigan un reparto ingenuo

Si el cubo raíz se corta en 8 octantes iguales, uno por trabajador, con 20.000 partículas queda así:

| Distribución | Partículas por octante | Máximo / promedio |
|---|---|---|
| uniforme | 2434, 2459, 2563, 2524, 2474, 2543, 2485, 2518 | 1,03 |
| dos_cumulos | 9866, 9950, 52, 52, 43, 32, 3, 2 | 3,98 |
| concentrada | 14219, 1691, 1430, 642, 583, 430, 414, 591 | 5,69 |

Todos los trabajadores esperan al más cargado (la barrera del final del paso), así que con dos cúmulos la aceleración máxima con 8 trabajadores sería 8 / 3,98 ≈ 2, no 8. Además, las partículas de un cúmulo denso abren más celdas del árbol, así que cada una cuesta más que una partícula aislada. Por eso el reparto debe hacerse por cantidad de trabajo (Morton, bisección recursiva o reparto dinámico) y no por volumen.

### La suma directa

Para cada partícula i:

```
a_i = Σ_j m_j (x_j − x_i) / (|x_j − x_i|² + ε²)^(3/2)
U   = −Σ_{i<j} m_i m_j / √(|x_j − x_i|² + ε²)
K   = Σ_i m_i |v_i|² / 2
```

Se calcula con numpy por **bloques de filas**: se toman unas cuantas partículas i y se calculan sus interacciones con todas las j a la vez. Si se hiciera de una sola vez, la matriz de diferencias tendría N × N × 3 números: para N = 20.000 serían 9,6 GB. Con bloques de 2 millones de pares, cada arreglo ocupa unos 50 MB sin importar N.

La autointeracción (i = j) se pone en cero a mano. Con ε > 0 ya daría cero, porque la diferencia es cero, pero con ε = 0 sería 0/0.

En la energía potencial, cada par se cuenta dos veces (desde i y desde j), por eso se multiplica por −1/2.

### El error de fuerzas

Por partícula: `|a_BH − a_directa| / |a_directa|`. Se reportan el promedio, la mediana y el máximo, junto con el índice de la peor partícula. El máximo sirve para encontrar un caso puntual que falla, por ejemplo una partícula en el borde de una celda.

### La simulación de referencia

`referencia_directa.py simular` integra con leapfrog usando fuerzas exactas. Sirve para separar dos fuentes de error en la energía:

- el error del **integrador**, que depende de Δt y existe aunque las fuerzas sean exactas;
- el error de **Barnes-Hut**, que depende de θ.

Si una versión tiene un error de energía mucho mayor que el de la referencia con el mismo Δt, el problema está en sus fuerzas.

Las aceleraciones del estado inicial se calculan **una vez antes del ciclo**, porque la primera media patada las necesita. Dentro del ciclo, las aceleraciones calculadas en el paso 3 se usan en la patada del paso 4 y también en la primera patada del paso siguiente. Así cada paso calcula fuerzas una sola vez.

## Resultados de las pruebas (en WSL)

| Prueba | Resultado |
|---|---|
| Dos partículas de masa 1 a distancia 1, ε = 0 | a = (±1, 0, 0) y U = −1, exacto |
| Energía potencial del cubo uniforme, N = 2.000, ε = 0,01 | −0,9422 (teórico para un cubo continuo: −0,9412) |
| Cálculo por bloques vs sin bloques, N = 3.000 | diferencia de 0 en aceleraciones y de 1e-9 en U |
| Verificar la referencia contra sí misma | error promedio 1,4e-10 (es el redondeo a 10 cifras del archivo) |
| Aceleraciones con ruido de 0,5 % / 5 % | aprobado / rechazado |
| Simulación directa, N = 2.000, Δt = 0,001, ε = 0,01, 100 pasos | error de energía 3,4e-8 (uniforme), 4,8e-5 (dos cúmulos), 5,9e-5 (concentrada) |
| Regenerar los archivos de 2.000 | idénticos (mismo md5) |

La suma directa tarda unos 0,3 s por cálculo de fuerzas con N = 2.000. Con 20.000 serían unos 30 s y con 200.000, casi una hora: por eso la verificación se hace con N pequeño.

## Decisiones de diseño

- **Una sola herramienta de verificación para las cuatro versiones.** Cada versión solo tiene que escribir dos archivos de texto. Así el criterio de correctitud es el mismo para todas y no puede haber una versión que se verifique "a su favor".
- **La energía se recalcula en `verificar.py`**, no se confía en la que imprime cada versión. Si una versión calcula mal la energía, igual se detecta.
- **Cúmulos en equilibrio y no en reposo.** Un cúmulo en reposo colapsa de golpe hacia el centro y la energía se dispara por el colapso mismo, no por un error del código. Con cúmulos en equilibrio, el error de energía refleja solo el integrador y la aproximación.
- **Tolerancias por defecto: 1 % en el error promedio de fuerzas y 0,1 % en la energía.** Con θ = 0,5 Barnes-Hut da errores promedio del orden de 0,5 % (medido en la fase 2), así que la tolerancia de 1 % está pensada para θ ≤ 0,5. La referencia exacta conserva la energía mejor que 1e-4. Pasar estos límites indica un error en el código. Las dos tolerancias se pueden cambiar con `--tolerancia`.
- **Diez cifras significativas en los archivos.** El redondeo (1e-10) queda muy por debajo de cualquier error que se mide, y los archivos son un 40 % más pequeños que con 17 cifras.
- **Masa total 1 en todas las distribuciones.** Así las energías de distintas distribuciones son comparables y las velocidades tienen la misma escala.

## Preguntas tipo sustentación

1. **¿Cómo saben que su versión de Barnes-Hut es correcta?**
   Escribe las aceleraciones del estado inicial sobre 2.000 partículas, y `verificar.py` las compara con la suma directa. Se reporta el error relativo promedio y máximo. Además, al final de la simulación se recalcula E = K + U de forma exacta y se compara con la inicial.

2. **¿Por qué se verifica con N = 2.000 y no con 200.000?**
   Porque la suma directa es O(N²). Con 2.000 partículas tarda 0,3 s y con 200.000 tardaría casi una hora por cada cálculo de fuerzas. El árbol no cambia de lógica con N, así que un N pequeño basta para encontrar errores.

3. **Si θ = 0, ¿qué error debería dar `verificar.py`?**
   Prácticamente cero, del orden del redondeo (1e-10 por las cifras del archivo). Con θ = 0 ninguna celda cumple el criterio, todas se abren hasta las hojas y el cálculo es una suma directa en otro orden.

4. **¿Por qué la energía no se conserva exactamente aunque las fuerzas sean exactas?**
   Porque leapfrog avanza en pasos finitos de Δt. Su error de energía no crece sin límite, sino que oscila, y se reduce como Δt². Con Δt = 0,001 queda entre 1e-8 y 1e-4.

5. **¿Qué es una esfera de Plummer y por qué se usa?**
   Un modelo de cúmulo con mucha densidad en el centro y poca afuera, que se sabe poner en equilibrio. Se usa para crear cúmulos densos que no colapsan de golpe, así la conservación de la energía mide el error del código.

6. **¿Por qué la distribución de dos cúmulos castiga un reparto ingenuo?**
   Si el cubo se divide en 8 partes iguales, dos partes tienen el 99 % de las partículas. Los demás trabajadores terminan casi enseguida y esperan en la barrera. Con 8 trabajadores, la aceleración no pasa de 2.

7. **¿Por qué la suma directa se calcula por bloques?**
   Porque la matriz completa de diferencias tiene N × N × 3 números (9,6 GB con N = 20.000). Por bloques se procesan unos 2 millones de pares a la vez y la memoria queda fija en unos 50 MB.

8. **¿Por qué se barajan las filas del archivo?**
   Para que ninguna versión dependa del orden del archivo. Si las partículas de un cúmulo vinieran juntas, repartir por índice daría un buen balanceo por casualidad. Barajadas, cada versión tiene que ordenar las partículas por posición.

9. **¿Por qué el leapfrog calcula aceleraciones antes del ciclo?**
   La primera media patada necesita las aceleraciones del estado inicial. Después, las aceleraciones del final de un paso sirven para la primera patada del siguiente, así cada paso calcula fuerzas una sola vez.

10. **Si el error promedio de fuerzas es bajo pero el máximo es muy alto, ¿qué puede pasar?**
    Un caso particular falla: por ejemplo, una partícula justo en el borde de una celda, una celda que se trata como lejana estando la partícula dentro de ella, o la partícula interactuando consigo misma. `verificar.py` imprime el índice de la peor partícula para revisarla.
