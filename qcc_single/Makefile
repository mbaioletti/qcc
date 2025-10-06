FLAGS=-O2

compiler: compiler.cpp problem.o solution.o
	g++ $(FLAGS) -o compile_qasm compiler.cpp problem.o solution.o
	
libcirc.so: circuito.cpp
	g++ -c -fPIC circuito.cpp -o circuito.o
	g++ -shared -Wl,-soname,libcirc.so -o libcirc.so circuito.o

problem.o: problem.cpp problem.h
	g++ $(FLAGS) -c problem.cpp
	
solution.o: solution.cpp solution.h problem.h
	g++ $(FLAGS) -c solution.cpp

clean:
	rm *.o compile_qasm
	
zip:
	zip qcc_v10.zip *.cpp *.h Makefile
