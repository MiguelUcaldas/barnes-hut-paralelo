# Fase 0: base del repositorio

## Qué se construyó

- La estructura de carpetas: una por versión (`python_mp/`, `mpi/`, `c_openmp/`, `cuda/`), una para lo compartido (`comun/`) y las de apoyo (`scripts/`, `datos/`, `resultados/`, `docs/`).
- `.gitattributes`: obliga a que todos los archivos de texto queden con saltos de línea LF.
- `.gitignore`: deja fuera del repositorio el entorno virtual, los compilados y los archivos grandes que se pueden regenerar.
- `requirements.txt`: numpy y matplotlib con versión fija.
- `README.md` y `ENTORNO.md`: qué es el proyecto, cómo preparar el entorno y con qué versiones se midió.

## Cómo funciona y por qué

**Saltos de línea.** Windows termina las líneas con CRLF y Linux con LF. Editamos en Windows y ejecutamos en WSL. Si un script de bash queda con CRLF, Linux lee el `\r` como parte del comando y falla con errores raros como `$'\r': command not found`. Con `* text=auto eol=lf`, git guarda y entrega todo con LF.

**Carpetas vacías.** Git no registra carpetas vacías, solo archivos. Por eso cada carpeta tiene un `.gitkeep` vacío mientras no tenga código.

**Versiones fijas.** El informe compara tiempos. Si cambia la versión de numpy, pueden cambiar los tiempos de Python, y entonces la medición no se podría repetir. Por eso `requirements.txt` usa `==`.

**Qué no se sube.** Los ejecutables (`bin/`, `*.o`) se generan al compilar y dependen de la máquina. Los archivos de 200.000 partículas pesan varios MB y se regeneran con semilla fija, así que se guardan en `datos/grandes/`, que está ignorada. En `datos/` solo quedan archivos pequeños de prueba.

**Cuatro versiones independientes.** Cada versión es un programa completo y no depende de las otras. Si una falla en la demostración, las demás siguen funcionando, y cada una se puede explicar por separado. Lo único compartido (en `comun/`) son las herramientas que no hacen parte de la simulación: generar las partículas, calcular la referencia exacta, verificar y graficar.

## Preguntas tipo sustentación

1. **¿Por qué el repositorio tiene `.gitattributes`?**
   Para que todos los archivos queden con LF. Editamos en Windows (CRLF) y ejecutamos en Linux. Un script de bash con CRLF falla en Linux.

2. **¿Por qué se fijan las versiones en `requirements.txt`?**
   Para que las mediciones se puedan repetir. Otra versión de numpy puede cambiar los tiempos de la versión de Python.

3. **¿Por qué los archivos grandes de partículas no se suben?**
   Porque se pueden regenerar idénticos con la semilla fija. Subirlos solo haría pesado el repositorio sin aportar información.

4. **¿Por qué las cuatro versiones no comparten el código de la simulación?**
   Porque el enunciado pide que cada versión resuelva su propia dificultad de paralelización. Además, si fueran independientes, un fallo en una no afecta a las demás.

5. **¿Qué significa que el equipo tenga 6 núcleos y 12 hilos?**
   Cada núcleo físico ejecuta dos hilos con SMT. Los dos hilos comparten las unidades de cálculo y la caché del núcleo. Por eso, pasar de 6 a 12 hilos rinde mucho menos que pasar de 3 a 6.

6. **¿Por qué `mpirun -np 8` falla por defecto en este equipo y cómo se arregla?**
   OpenMPI da un slot por núcleo físico, o sea 6. Con `--use-hwthread-cpus` cuenta los 12 hilos. Con `--oversubscribe` permite más procesos que slots.

7. **¿Qué se espera al usar 16 trabajadores en un equipo de 12 hilos?**
   No mejora y suele empeorar. Hay más trabajadores que hilos de hardware, así que el sistema operativo los turna, y eso agrega cambios de contexto y competencia por la caché.

8. **¿Por qué WSL2 ve solo unos 7,4 GB de los 16 GB?**
   WSL2 corre en una máquina virtual que por defecto recibe la mitad de la memoria del equipo. Esto limita el N máximo que se puede simular, sobre todo en MPI, donde cada proceso tiene su propia copia de los datos.
