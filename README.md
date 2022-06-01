# Mesh-Fairing
A C++ program to read a file in a simplified obj format and write an obj file of a smoothed mesh.

How to run:-

First run the make make command to smooth
then use any flags to generate desired image

Basic Assignment:-
 
Use:	```make smooth ```
	
to generate out file 
the output is "smooth"

Use:	```smooth bunny.obj output.obj 1 50```
 	for basic output, here 1 is the step size and 50 is the iterations, 
	and bunny obj is input file. 

 the output is "output.obj"

For Implicit Integration:-
Use:	```smooth bunny.obj output.obj 1 50 -i 0.000005 ```                  [must specify eps only at the end]
 	for basic output, here 1 is the step size and 50 is the iterations, and -i is to specify eps 
	and bunny obj is input file. 

 the output is "output.obj"
      

Note:-
Implimented Extra Credit Implicit Itegration, please make sure you specify eps at the end with -i, not in the middle or at the front
 
Resources used:-

Program Webiste 
https://www.csee.umbc.edu/~adamb/435/proj4.html
https://www.csee.umbc.edu/~adamb/435/implicit.cpp // for extra credit Implict integration

Paper the Assig was based on:-
http://multires.caltech.edu/pubs/ImplicitFairing.pdf

LinearSolver from Glserver for Extra Credit Implicit Integration
