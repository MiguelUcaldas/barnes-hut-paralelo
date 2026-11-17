# Compilar desde esta carpeta: make
# -std=gnu17 es C17 más las funciones POSIX que se usan (clock_gettime, getopt_long, sysconf).
CC = gcc
CFLAGS = -O2 -std=gnu17 -Wall -Wextra -Wpedantic
LDLIBS = -lm

FUENTES = barnes_hut.c configuracion.c particulas.c arbol.c fuerzas.c medicion.c
CABECERAS = configuracion.h particulas.h arbol.h fuerzas.h medicion.h memoria.h

bin/barnes_hut: $(FUENTES) $(CABECERAS)
	mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $(FUENTES) $(LDLIBS)

# Detecta accesos fuera de los arreglos, memoria liberada y comportamiento
# indefinido. Es varias veces más lenta: solo para probar con N pequeño.
bin/barnes_hut_sanitizado: $(FUENTES) $(CABECERAS)
	mkdir -p bin
	$(CC) -O1 -g -std=gnu17 -Wall -Wextra -fsanitize=address,undefined -o $@ $(FUENTES) $(LDLIBS)

sanitizado: bin/barnes_hut_sanitizado

clean:
	rm -rf bin

.PHONY: sanitizado clean
