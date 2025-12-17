FLAGS=-O2

compile_qasm: main.cpp compiler.o problem.o solution.o
	g++ $(FLAGS) -o compile_qasm main.cpp compiler.o problem.o solution.o
	
compiler.o: compiler.cpp compiler.h
	g++ $(FLAGS) -c compiler.cpp 
	
libcirc.so: compile_from_qiskit.cpp
	g++ -c -fPIC compile_from_qiskit.cpp -o compile_from_qiskit.o
	g++ -c -fPIC compiler.cpp -o compiler_qsk.o
	g++ -c -fPIC problem.cpp -o problem_qsk.o
	g++ -c -fPIC solution.cpp -o solution_qsk.o
	g++ -shared -Wl,-soname,libcirc.so -o libcirc.so compile_from_qiskit.o problem_qsk.o solution_qsk.o compiler_qsk.o

problem.o: problem.cpp problem.h
	g++ $(FLAGS) -c problem.cpp
	
solution.o: solution.cpp solution.h problem.h
	g++ $(FLAGS) -c solution.cpp

clean:
	rm *.o compile_qasm
	
zip:
	zip qcc_d_3.zip *.cpp *.h Makefile
