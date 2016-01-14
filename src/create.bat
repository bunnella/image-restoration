#!/usr/bin/bash
make -C lib/png
gcc -o test test.c -Llib/png -lpng -lz
