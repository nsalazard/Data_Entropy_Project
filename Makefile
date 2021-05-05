# Esto es un comentario
# @: target
# ^: todas las dependencias
# <: la primera dependencia
CXXFLAGS = -lpapi -larmadillo
CXXO = -O3

all: datos.txt

%.x : %.cpp
	g++ $< -o $@ $(CXXFLAGS)


#CPU vs Nb

#Sin optimizar
mult_2048.txt : multiplicacion_blocking.x
	./$< > $@ 0 0

#Con optimización
%_O.x : %.cpp
	g++ $(CXXO) $< -o $@ $(CXXFLAGS)

mult_2048.txt : multiplicacion_blocking_0.x
	./$< > $@ 0 0





#Crecaion general de graficas
%.pdf: %.gp %.txt 
	gnuplot $<
	
clean:
	rm -f a.out *.x *.txt
