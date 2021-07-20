#staticLink = ./lib/libcblas_WIN.a ./lib/*.o
include = -I ./include -L ./lib
remove = rm
opts = -fopenmp
options = ${include} $(opts) #-Wall

argsO = "-O3" "-g" -D OMP_ON
argsd = "-g"

main.exe: main.cpp *.hpp *.h
	g++ main.cpp -o main.exe -g -fopenmp ${argsd} ${options}

testRoe.exe: testRoe.cpp
	g++ testRoe.cpp -g -o testRoe.exe ${argsd} ${options}

mainSG.exe: mainSG.cpp roesgsolver.hpp
	g++ mainSG.cpp -g -o mainSG.exe ${argsd} ${options}

release: mainR.exe

mainR.exe: main.cpp *.hpp *.h
	g++ main.cpp -o mainR.exe -g -fopenmp ${argsO} ${options}

mainSGR.exe: mainSG.cpp roesgsolver.hpp
	g++ mainSG.cpp -g -o mainSGR.exe ${argsO} ${options}

.PHONY: clean
clean:
	$(remove) *.exe   