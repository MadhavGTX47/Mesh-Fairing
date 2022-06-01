
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "argument.h"

#include "lsolver/cghs.h"
#include "lsolver/bicgsq.h"
#include "lsolver/bicgstab.h"
#include "lsolver/gmres.h"

using std::cout;
using std::endl;
using std::cerr;


// Steifigkeitsmatrix
struct StiffMat {
  StiffMat( int N, double h ) : N(N), h(h) { }
  int N;
  double h;
};

// w=A*v  (Steifigkeitsmatrix A)
void mult( const StiffMat &A, const double *v, double *w ) {
  int n=A.N;
  int Noff = n*(n+1);
  for ( int off=n+1; off<Noff; off+=n+1 ) {
    int Ni = off+n;
    for ( int i=off+1; i<Ni; ++i )
      w[i] = (4.*v[i]-v[i-n-1]-v[i+n+1]-v[i-1]-v[i+1]);
  }
}


// Massenmatrix
struct MassMat {
  MassMat( int N, double h ) : N(N), h(h) { }
  int N;
  double h;
};

// w=M*v  (Massenmatrix M)
void mult( const MassMat &M, const double *v, double *w ) {
  double h=M.h;
  int n=M.N;
  int Noff = n*(n+1);
  for ( int off=n+1; off<Noff; off+=n+1 ) {
    int Ni = off+n;
    for ( int i=off+1; i<Ni; ++i )
      w[i] = h*h*(.5*v[i]+(v[i-n-1]+v[i+n+1]+v[i-1]+v[i+1]+v[i+n+2]+v[i-n-2])/12.);
  }
}



//
// Aufgabe: -laplace u=f mit L"osung g
//

const double pi = 3.14159265358979323846;


#if 1
inline
double f( double x, double y ) {
  return 8.*pi*pi*sin(2.*pi*x)*sin(2.*pi*y);
}

inline
double g( double x, double y ) {
  return sin(2.*pi*x)*sin(2.*pi*y);
}
#else
inline
double f( double x, double y ) {
  return pi*pi*sin(pi*x)*sin(pi*y);
}

inline
double g( double x, double y ) {
  return sin(pi*x)*sin(pi*y);
}
#endif





// call: serial N   =>  discretisation with {0,1,...,N}x{0,1,...,N}
int main( int argc, char *argv[] ) {
  Arguments args(argc,argv);
  int N         = args.get_int_argument("-n",50);  // #Punkte in eine Richtung
  int m         = args.get_int_argument("-m",5);   // f"ur gmres(M)
  double EPS    = args.get_double_argument("-eps",1e-5);
  bool CGHS     = args.contains_flag("-cghs");
  bool BICGSQ   = args.contains_flag("-bicgsq");
  bool BICGSTAB = args.contains_flag("-bicgstab");
  bool GMRES    = args.contains_flag("-gmres");
  if ( !args.no_more_args() ) {
    cerr<<"falsche Argumente"<<endl;
    exit(-1);
  }

  int start = clock() ;

  double h=1./(double)N;         // grid width
  int M=(N+1)*(N+1);             // #nodes

  double *u = new double[M];
  double *v = new double[M];
  double *w = new double[M];
  double *b = new double[M];
  {
    int i,j,k;
    double x,y;
    for ( k=0; k<M; ++k )
      u[k]=v[k]=w[k]=b[k]=0.;
    k=N+2;
    for ( i=1,x=h; i<N; ++i,x+=h ) {
      for ( j=1,y=h; j<N; ++j, y+=h ) {
	// k = (N+1)*i+j
	u[k]=0.;       // L"osungsvektor
	v[k]=g(x,y);   // exakte L"osung
	w[k]=f(x,y);   // -laplace(u)=f
	++k;
      }
      k+=2;
    }
  }
  StiffMat  A(N,h);
  MassMat   B(N,h);
  mult(B,w,b);                 // berechne b f"ur Au=b (A Steifigkeitsmatrix)

  int its = -1;
  if ( CGHS )
    its = cghs(M,A,b,u,EPS);
  else if ( BICGSQ )
    its = bicgsq(M,A,b,u,EPS);
  else if ( BICGSTAB )
    its = bicgstab(M,A,b,u,EPS);
  else if ( GMRES )
    its = gmres(m,M,A,b,u,EPS);
  else {
    cerr<<"kein L\"oser (-cghs -bicgsq -bicgstab -gmres) angegeben"<<endl;
    exit(-1);
  }

  // Differenz zur exakten L"osung
  {
    for ( int i=0; i<M; ++i )
      w[i] = u[i]-v[i];
  }
  // maximale Abweichung
  double max = 0;
  {
    for ( int i=0; i<M; ++i )
      max= ( fabs(w[i])>max ) ? fabs(w[i]) : max;
  }
  cout<<its<<" iterations, difference="<<max<<endl;

  float used = (float)(clock() - start)/(float)(CLOCKS_PER_SEC);
  cout<<"needed time: "<<used<<endl;

  delete[] u;
  delete[] v;
  delete[] w;
  delete[] b;
}
