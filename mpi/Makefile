# Compilar: make -C mpi (desde la raíz) o make (desde esta carpeta).
# mpicc es gcc con las rutas de las cabeceras y bibliotecas de OpenMPI.
# -std=gnu17 es C17 más las funciones POSIX que se usan (getopt_long, sysconf).
CC = mpicc
CFLAGS = -O2 -std=gnu17 -Wall -Wextra
LDLIBS = -lm

FUENTES = barnes_hut.c configuracion.c entrada_salida.c particion.c morton.c arbol.c esencial.c fuerzas.c medicion.c
CABECERAS = configuracion.h entrada_salida.h particion.h morton.h arbol.h criterio.h esencial.h fuerzas.h medicion.h particulas.h memoria.h

bin/barnes_hut: $(FUENTES) $(CABECERAS) Makefile
	mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $(FUENTES) $(LDLIBS)

# Detecta comportamiento indefinido (desbordes de enteros, desplazamientos de
# bits inválidos, accesos desalineados) y se detiene en el primero. No se usa
# AddressSanitizer: en este sistema choca con OpenMPI dentro de MPI_Init,
# incluso con un programa MPI vacío. Más lenta: solo para probar con N pequeño.
bin/barnes_hut_sanitizado: $(FUENTES) $(CABECERAS) Makefile
	mkdir -p bin
	$(CC) -O1 -g -std=gnu17 -Wall -Wextra -fsanitize=undefined -fno-sanitize-recover=all -o $@ $(FUENTES) $(LDLIBS)

sanitizado: bin/barnes_hut_sanitizado

clean:
	rm -rf bin

.PHONY: sanitizado clean
