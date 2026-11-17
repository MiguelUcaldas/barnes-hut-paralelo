# Compilar: make -C c_openmp (desde la raíz) o make (desde esta carpeta).
# -std=gnu17 es C17 más las funciones POSIX que se usan (clock_gettime, getopt_long, sysconf).
# -fopenmp activa los #pragma omp y enlaza la biblioteca de OpenMP.
CC = gcc
CFLAGS = -O2 -std=gnu17 -fopenmp -Wall -Wextra -Wpedantic
LDLIBS = -lm

FUENTES = barnes_hut.c configuracion.c particulas.c morton.c arbol.c fuerzas.c medicion.c
CABECERAS = configuracion.h particulas.h morton.h arbol.h fuerzas.h medicion.h memoria.h

bin/barnes_hut: $(FUENTES) $(CABECERAS) Makefile
	mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $(FUENTES) $(LDLIBS)

# Detecta accesos fuera de los arreglos, memoria liberada y comportamiento
# indefinido. Es varias veces más lenta: solo para probar con N pequeño.
bin/barnes_hut_sanitizado: $(FUENTES) $(CABECERAS) Makefile
	mkdir -p bin
	$(CC) -O1 -g -std=gnu17 -fopenmp -Wall -Wextra -fsanitize=address,undefined -o $@ $(FUENTES) $(LDLIBS)

sanitizado: bin/barnes_hut_sanitizado

clean:
	rm -rf bin

.PHONY: sanitizado clean
