CFLAGS=-std=c99 -O3
CC=gcc $(CFLAGS)

.PHONY: clean test

all: runEMD

runEMD: emd.o runEMD.o
	$(CC) runEMD.o emd.o -lm -o runEMD

emd.o: emd.c emd.h
	$(CC) -c emd.c -o emd.o

runEMD.o: runEMD.c emd.h
	$(CC) -c runEMD.c -o runEMD.o

clean:
	rm -v *.o runEMD

test:
	mkdir -pv test/run
	./runEMD test/1L5BA_sig.txt test/1L5EA_sig.txt 288 > test/run/1L5BA-1L5EA.out
	diff test/expected/1L5BA-1L5EA.out test/run/1L5BA-1L5EA.out
