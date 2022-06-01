smooth : io.o smooth.o 
	g++ smooth.o io.o -lblas -o smooth

smooth.o : smooth.cpp 
	g++ -I/usr/include/eigen3/ -O3 -c smooth.cpp
io.o :
	g++ -I/usr/include/eigen3/ -O3 -c io.cpp

clean:
	rm *.o *.ppm *.exe output out
	 
