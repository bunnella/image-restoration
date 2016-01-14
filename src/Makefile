all: libpng test

libpng: 
	make -C lib/png

test: 
	gcc -o test test.c -Llib/png -lpng -lz
