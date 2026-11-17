# Barnes-Hut paralelo

Simulación gravitacional de N cuerpos en 3D (20.000 a 200.000 partículas) con el algoritmo de Barnes-Hut, paralelizada en cuatro versiones independientes:

| Carpeta | Versión | Cómo se paraleliza |
|---|---|---|
| `python_mp/` | Python con `multiprocessing` | Árbol en arreglos planos de numpy, compartido o reconstruido entre procesos |
| `mpi/` | C con MPI (OpenMPI) | Particionado del espacio y envío solo de las partes remotas del árbol que cada proceso necesita |
| `c_openmp/` | C con OpenMP | Construcción paralela del árbol y balanceo de carga entre hilos |
| `cuda/` | CUDA en Google Colab | Árbol aplanado, partículas ordenadas y recorrido sin recursión |

Cada versión es un programa completo: construye el árbol, calcula las fuerzas e integra. Lo compartido son las herramientas de `comun/` (generar partículas, suma directa de referencia, verificación y gráficas).

Proyecto final de Programación Concurrente y Distribuida, Universidad de Caldas. El enunciado está en [docs/enunciado-proyecto-final.pdf](docs/enunciado-proyecto-final.pdf).

## Física y convenciones

- G = 1 y suavizado ε: la aceleración usa (r² + ε²)^(3/2) en el denominador y la energía potencial usa √(r² + ε²).
- Integrador leapfrog patada-deriva-patada: medio paso de velocidad, paso completo de posición, recálculo de aceleraciones con Barnes-Hut y otro medio paso de velocidad.
- Criterio de apertura: una celda se trata como una sola partícula si `tamaño_celda / distancia < θ`; si no, se abre. Con θ = 0 el resultado es igual a la suma directa.
- El árbol se reconstruye en cada paso a partir del cubo más pequeño que contiene todas las partículas.

## Formato de las partículas

Texto plano, una partícula por fila, siete columnas separadas por espacios:

```
x y z vx vy vz masa
```

## Estructura

```
comun/        herramientas compartidas en Python
c_openmp/     versión C + OpenMP
mpi/          versión MPI
python_mp/    versión Python + multiprocessing
cuda/         notebook para Google Colab
scripts/      medición automática
datos/        archivos de partículas de prueba
resultados/   mediciones (CSV) y gráficas
docs/         enunciado, mapa del código y explicaciones por fase
```

## Preparar el entorno

Se ejecuta en Ubuntu (WSL2). Los detalles y las versiones exactas están en [ENTORNO.md](ENTORNO.md).

```bash
sudo apt install build-essential openmpi-bin libopenmpi-dev htop python3-venv
git clone https://github.com/MiguelUcaldas/barnes-hut-paralelo.git
cd barnes-hut-paralelo
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Generar las partículas

```bash
python comun/generar_particulas.py --distribucion uniforme --n 2000      # uniforme, dos_cumulos o concentrada
bash scripts/generar_datos.sh .venv/bin/python                           # todos: 2.000 en datos/, 20.000 a 200.000 en datos/grandes/
```

Los archivos de 2.000 partículas vienen en `datos/`. Los grandes no se suben a git; se regeneran idénticos con la semilla fija.

## Verificar una versión

Cada versión acepta `--guardar-aceleraciones archivo` (aceleraciones del estado inicial, `ax ay az` por fila) y `--guardar-final archivo` (estado final, 7 columnas). Con esos archivos:

```bash
python comun/verificar.py fuerzas --particulas datos/uniforme_2000.txt --aceleraciones acel.txt --epsilon 0.01
python comun/verificar.py energia --inicial datos/uniforme_2000.txt --final final.txt --epsilon 0.01
```

La primera compara contra la suma directa O(N²) y reporta el error relativo promedio, mediano y máximo. La segunda recalcula E = K + U de forma exacta al inicio y al final. Para conocer el error de energía que pone el integrador solo, sin Barnes-Hut:

```bash
python comun/referencia_directa.py simular --particulas datos/uniforme_2000.txt --epsilon 0.01 --dt 0.001 --pasos 100
```

## Compilar y ejecutar

Todos los programas se ejecutan desde la raíz del repositorio.

### C (`c_openmp/`)

```bash
make -C c_openmp                 # bin/barnes_hut
make -C c_openmp sanitizado      # bin/barnes_hut_sanitizado: detecta errores de memoria (lento)
c_openmp/bin/barnes_hut --entrada datos/uniforme_2000.txt --pasos 100 --theta 0.5 --epsilon 0.01 --dt 0.001
c_openmp/bin/barnes_hut --ayuda
```

Opciones: `--entrada`, `--pasos`, `--dt`, `--theta`, `--epsilon`, `--cada-energia K`, `--guardar-aceleraciones f`, `--guardar-final f`, `--distribucion nombre`, `--csv f`.

## Resultados

Cada versión imprime al terminar la verificación (error de fuerzas contra la suma directa y conservación de la energía) y los tiempos por fase. También agrega una fila a `resultados/mediciones.csv` con estas columnas:

```
version, distribucion, N, pasos, theta, epsilon, dt, trabajadores, t_construccion, t_ordenamiento, t_fuerzas, t_integracion, t_total, energia_inicial, energia_final, error_energia
```
