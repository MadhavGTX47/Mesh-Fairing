#include <string>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <Eigen/Geometry>
#include <assert.h>
#include "io.h"    
#include "LinearSolver/cghs.h"

using namespace Eigen;
using namespace std;

struct Operator {  // Code from implicit.cpp from course site
  int n;
  double  dt;
  vector<Tri> tris;
};


int smoothit(char* input_obj, char* output_obj, double stepsize, double iteration);	// Declaring smooth function	
int implicit_integrateit(char* input_obj, char* output_obj, double stepsize, double eps); // Declaring Extra credit Implicit integration
void mult (const Operator &op, double *v, double *w);
vector<Vector3d> laplacian(const vector<Vector3d> &pts, const std::vector<Eigen::Vector3d> &tris);





char *input_obj;	//our input and output
char *output_obj;
double eps=-1; 

int main(int argc, char* argv[])
{
	if (argc < 4) {           
        cout<<"Check number of Arguments"; // checking the number of arguments
        exit(1);
    }

	input_obj = argv[1];	//Storing out Arguments
	output_obj = argv[2];
	double stepsize = stod(argv[3]);
	double iteration = stod(argv[4]);

	int c;
	 while (( c = getopt (argc, argv, "i:")) != -1)  //Handling command line flags
    switch (c)
      {
      case 'i':
       eps=stod(optarg);
        break;
      
      default:
        abort ();
      }

	if(eps==-1)
    smoothit(input_obj,output_obj,stepsize,iteration); //Calling smooth function
	
	else
	implicit_integrateit(input_obj,output_obj,stepsize,eps); //Calling Implicit Integration function

return 0;
}		


int smoothit(char* input_obj, char* output_obj, double stepsize, double iteration){ //Smoothing base assig function

	vector<Vector3d> vertex;  //Decalaring variables for traingels and vertices
	vector<Tri> triangles;
	vector<double> m;
	vector<Eigen::Vector3d> L;
	readObjFile(input_obj, vertex, triangles);

	for(int k=0; k< iteration;k++){ //Main loop for number of iterations

		for (unsigned int j=0; j<vertex.size(); j++) { // looping through vertices
			
			Vector3d temp= {0.0,0.0,0.0}; //setting L and m to 0
            L.push_back(temp);
            m.push_back(0);
     	}

		  for (unsigned int i=0; i<triangles.size(); i++) { //Umbrella Operation and looping through triangles
            Tri &tri = triangles[i];      //storing triangle at i to tri
			
            for (int z=0; z<3; z++) {
                L[tri.indices[z]] += vertex[tri.indices[(z+1)%3]] + vertex[tri.indices[(z+2)%3]] - 2 * vertex[tri.indices[z]]; //calculating L[tri[0]]L[tri[1]]L2[tri[3]]]
				m[tri.indices[z]] += 2; 
            }

        }

		for (unsigned int i=0; i<vertex.size(); i++) { //looping over vertices
            vertex[i] +=  stepsize * (L[i]/m[i])  ;    //Updating vertices positions
        }

	}

		
    writeObjFile(output_obj, vertex, triangles); //Writing output to Object file

    cout<<"done";
    return 0;


}

int implicit_integrateit(char* input_obj, char* output_obj, double stepsize, double eps){ // Used Sample code from implicit.cpp from course site
	
	struct Operator op; //Decalaring variables for vertices and Op
	vector<Vector3d> vertex;
  
	
	readObjFile(input_obj, vertex, op.tris); // Reading the the input obj file

	op.n = vertex.size();  op.dt = stepsize;

	double b[3*vertex.size()];  //Initializing b and x
  	double x[3*vertex.size()];



	for (unsigned int i=0; i<vertex.size(); i++) {   //Assiging values for b and x from vertices
        for (int c=0; c<3; c++) {
            b[3*i+c] = vertex[i][c];
            x[3*i+c] = vertex[i][c];
        }
    }
	
	cghs<Operator>(3*vertex.size(), op, b, x, eps, true); // Calling the cghs in LinerSolver

	for (unsigned int i=0; i<vertex.size(); i++) {
        vertex[i] = Vector3d(x[3*i+0], x[3*i+1], x[3*i+2]); // writing vertex with returned L
    }

    writeObjFile(output_obj, vertex, op.tris); //Writing output to Object file
    cout<<"done";
    return 0;
}


vector<Vector3d> laplacian(const vector<Vector3d> &vertex, const vector<Tri> &triangles) {
    std::vector<Eigen::Vector3d> L;
	
	for (unsigned int j=0; j<vertex.size(); j++) { // looping through vertices

			Vector3d temp= {0.0,0.0,0.0}; //setting L  to 0
            L.push_back(temp);
     	}

	 for (unsigned int i=0; i<triangles.size(); i++) { // looping through triangles
           const Tri &tri = triangles[i];      //storing triangle at i to tri
			
            for (int z=0; z<3; z++) {
                L[tri.indices[z]] += vertex[tri.indices[(z+1)%3]] + vertex[tri.indices[(z+2)%3]] - 2 * vertex[tri.indices[z]]; //calculating L[tri[0]],L[tri[1]],L2[tri[3]]] 
            }

        }

    return L; 
}

void mult (const Operator &op, double *v, double *w) { // Code from implicit.cpp from course site
  vector<Eigen::Vector3d> vertex;
  vertex.resize(op.n);
  for (unsigned int i=0; i<vertex.size(); i++) {
      vertex[i] = Eigen::Vector3d(v[3*i+0], v[3*i+1], v[3*i+2]);
  }

 vector<Vector3d> l = laplacian(vertex, op.tris);
  for (unsigned int i=0; i<vertex.size(); i++) {
    l[i] *= op.dt;
    w[3*i+0] = v[3*i+0] - l[i][0];
    w[3*i+1] = v[3*i+1] - l[i][1];
    w[3*i+2] = v[3*i+2] - l[i][2];
  }
}


		



