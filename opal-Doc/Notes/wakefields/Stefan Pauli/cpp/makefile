
all: run

run: wake
	./a.out

wake: wake.o twoParticle.o task.o
	g++ wake.o twoParticle.o task.o -o a.out  -O2  -march=pentium4 -fopenmp  -lfftw3 -lm


wake.o: wake.cpp
	g++ -c wake.cpp -O2  -march=pentium4 -fopenmp -lfftw3 -lm

twoParticle.o: twoParticle.cpp
	g++ -c twoParticle.cpp -O2  -march=pentium4 -fopenmp -lfftw3 -lm
	
task.o: task.cpp
	g++ -c task.cpp -O2  -march=pentium4 -fopenmp 

debug: wake
	./a.out

clear:
	rm wake.o -f
	rm twoParticle.o -f
	rm task.o -f

.phony: clear all

