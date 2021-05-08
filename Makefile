# Esto es un comentario
# @: target
# ^: todas las dependencias
# <: la primera dependencia
CXXFLAGS = -lpapi -larmadillo
CXXO = -O3

all :plot.pdf plot1.pdf plot2.pdf plot3.pdf

plot.pdf: plot.gp bloc_2.txt bloc_4.txt bloc_8.txt 
	gnuplot plot.gp

plot1.pdf: plot1.gp bloc_16.txt bloc_32.txt bloc_64.txt 
	gnuplot plot1.gp

plot2.pdf: plot2.gp bloc_128.txt bloc_256.txt 
	gnuplot plot2.gp 

plot3.pdf: plot3.gp bloc_512.txt bloc_1024.txt
	gnuplot plot3.gp

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
