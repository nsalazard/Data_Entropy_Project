# Esto es un comentario
# @: target
# ^: todas las dependencias
# <: la primera dependencia
CXXFLAGS = -lpapi -larmadillo
CXXO = -O3

all: mult_2048_O.txt

%.x : %.cpp
	g++ $< -o $@ $(CXXFLAGS)


#CPU vs Nb

#Sin optimizar
mult_2048.txt : multiplicacion_blocking.x
	./$< > $@ 0 0

#Con optimización
multiplicacion_blocking_O.x : multiplicacion_blocking.cpp
	g++ $(CXXO) $< -o $@ $(CXXFLAGS)

mult_2048_O.txt : multiplicacion_blocking_O.x
	./$< > $@ 0 0

	
clean:
	rm -f a.out *.x *.txt
