# Esto es un comentario
# @: target
# ^: todas las dependencias
# <: la primera dependencia
CXXFLAGS = -lpapi -larmadillo
CXXO = -O3

all: mult_2048_O.txt mult_2048.txt

%.x : %.cpp
	g++ $< -o $@ $(CXXFLAGS)


#CPU vs Nb

#Sin optimizar
mult_2048.txt : n.x
	./$< 0 0 > $@ 

#Con optimización
n_O.x : n.cpp
	g++ $(CXXO) $< -o $@ $(CXXFLAGS)

mult_2048_O.txt : n_O.x
	./$< 0 0 > $@

	
clean:
	rm -f a.out *.x *.txt
