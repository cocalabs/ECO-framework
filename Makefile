#!/bin/sh

#To use GDB put '-g' option in 'eco-cluster:'
#Run the executable as gdb ./name
#In gdb execute 'r input_eco.in' and than 'tb'

CC=gcc

#OBJECTS = example.o cluster.o
OBJECTS = mersenne.o cluster.o eco-cluster.o
SRCDIR = ./src-clstr

eco-cluster: $(OBJECTS)
	$(CC) $(OBJECTS) -o eco-cluster -Wall -Wextra -Wformat -Wformat-security -O3 -lpthread -lm

eco-cluster.o: eco-cluster.c eco.h
	$(CC) eco-cluster.c -c -O3 -lpthread -lm

cluster.o: $(SRCDIR)/cluster.c $(SRCDIR)/cluster.h
	$(CC) $(SRCDIR)/cluster.c -c -Wall -lm

mersenne.o: mersenne.c
	$(CC) mersenne.c -c -lm

clean:
	rm eco-cluster eco-cluster.o cluster.o mersenne.o

