#!/bin/bash
gcc spectralclustering_input_matrix_v2.c -o spectralclustering.o /usr/lib/liblapack.so -llapack -lm
