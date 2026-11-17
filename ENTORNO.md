# Entorno

## Equipo de la demostración

| | |
|---|---|
| Procesador | AMD Ryzen 5 5500U: 6 núcleos físicos, 12 hilos (SMT) |
| Memoria | 16 GB DDR4 en doble canal (WSL2 ve unos 7,4 GB) |
| Disco | SSD NVMe |
| Sistema anfitrión | Windows 11 |
| Sistema de ejecución | Ubuntu 26.04.1 LTS en WSL2 (kernel 6.18.33.1-microsoft-standard-WSL2) |

El código se edita en Windows y se compila y ejecuta dentro de WSL2.

## Versiones

| Herramienta | Versión |
|---|---|
| gcc (con OpenMP) | 15.2.0 |
| OpenMPI | 5.0.10 (PMIx 5.0.9, PRRTE 3.0.13) |
| Python | 3.14.4 |
| numpy | 2.5.5 |
| matplotlib | 3.11.2 |
| htop | 3.4.1 |
| GPU (Colab) | NVIDIA T4 |

## Replicar el entorno

```bash
sudo apt update
sudo apt install build-essential openmpi-bin libopenmpi-dev htop python3-venv
git clone https://github.com/MiguelUcaldas/barnes-hut-paralelo.git ~/barnes-hut-paralelo
cd ~/barnes-hut-paralelo
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

## Notas sobre MPI en este equipo

- Por defecto `mpirun` da un slot por núcleo físico, así que permite hasta 6 procesos.
- Para 8 a 12 procesos hay que contar los hilos SMT: `mpirun --use-hwthread-cpus -np 8 ...`
- Para más de 12 procesos: `mpirun --oversubscribe -np 16 ...`
- En Ubuntu 26.04 faltan los textos de ayuda de PRRTE, entonces los errores salen solo con su nombre. Por ejemplo, `prte-rmaps-base:alloc-error` significa que no hay suficientes slots para los procesos pedidos (se resuelve con las opciones de arriba).
