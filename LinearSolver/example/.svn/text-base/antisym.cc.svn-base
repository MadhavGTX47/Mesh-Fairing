
// 
// gegeben f(n), n=1,2,...
// Matrix A wird so erstellt, dass A[i,i]=1 und
//  A[i,j]=f[i+j] (i>j) sowie A[i,j]=-f[i+j] (i<j)
//



#include <iostream.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "argument.h"

#include "lsolver/cghs.h"
#include "lsolver/bicgsq.h"
#include "lsolver/bicgstab.h"
#include "lsolver/gmres.h"


//
// verschiedene Funktionen f
//
inline double pow( double a, int n ) {
  double x=1.;
  while ( n-- )
    x*=a;
  return x;
}

double vandermonde( int n ) {
  return 1./double(n);
}

double quadratic( int n ) {
  return 1./double(n*n);
}

double a=1.2;
double geom( int n ) {
  return 1./double(pow(a,n));
}




//
// erstellen der Matrix
//
struct AntiMat {
  double *x;     // Wert der k-ten Hauptdiagonalen
  double reyn;   // St"arke der Antisymmetrie
  int n;         // Dimension
};

void setupAntiMat( int n, double(*f)(int), double reyn, AntiMat &A ) {
  A.n = n;
  A.reyn = reyn;
  A.x = new double[n];
  for ( int i=1; i<n; ++i )
    A.x[i] = f(i);
}

void destructAntiMat( AntiMat &A ) {
  delete[] A.x;
}

//
// Multiplikation Matrix*Vektor
//
void mult( const AntiMat &A, const double *v, double *w ) {
  int n=A.n;
  // 1. die Diagonale
  for ( int k=0; k<n; ++k )
    w[k] = v[k];
  // untere/obere Dreiecksmatrix
  for ( int i=0; i<n; ++i ) {
    for ( int j=0; j<i; ++j ) {
      w[i] += A.reyn*A.x[i-j]*v[j];
      w[j] -= A.reyn*A.x[i-j]*v[i];
    }
  }
}




int main( int argc, char *argv[] ) {
  Arguments args(argc,argv);
  int N         = args.get_int_argument("-n",50);       // Dimension
  int M         = args.get_int_argument("-m",5);        // f"ur gmres(M)
  double EPS    = args.get_double_argument("-eps",1e-5);// iteriere bis Res<eps
  bool CGHS     = args.contains_flag("-cghs");
  bool BICGSQ   = args.contains_flag("-bicgsq");
  bool BICGSTAB = args.contains_flag("-bicgstab");
  bool GMRES    = args.contains_flag("-gmres");
  bool LIN      = args.contains_flag("-lin");           // Fkt vandermonde()
  bool QUAD     = args.contains_flag("-quad");          // Fkt quadratic()
  bool GEOM     = args.contains_flag("-geom");          // Fkt geom()
  a             = args.get_double_argument("-a",1.2);   // f"ur geom-Funktion
  double REYN   = args.get_double_argument("-reyn",1.); // St"arke d. Antisym.
  if ( !args.no_more_args() ) {
    cerr<<"falsche Argumente"<<endl;
    exit(-1);
  }

  int start = clock() ;
  double *u = new double[N];
  double *v = new double[N];
  double *w = new double[N];
  double *b = new double[N];

  { // rechte Seite = konstant 1
    for ( int i=0; i<N; ++i )
      b[i] = 1.;
  }
  { // Iterationsvektor mit 0 initialisieren
    for ( int i=0; i<N; ++i )
      u[i] = 0.;
  }

  AntiMat     A;
  if ( LIN )
    setupAntiMat(N,vandermonde,REYN,A);
  else if ( QUAD )
    setupAntiMat(N,quadratic,REYN,A);
  else if ( GEOM )
    setupAntiMat(N,geom,REYN,A);
  else {
    cerr<<"keine Belegungsfunktion (-lin -quad -geom) angegeben"<<endl;
    exit(-1);
  }

  // Gleichungssystem l"osen
  int its = -1;
  if ( CGHS )
    cout<<"Achtung: CGHS f\"ur nicht-symmetrische Werte nicht definiert"<<endl,
    its = cghs(N,A,b,u,EPS);
  else if ( BICGSQ )
    its = bicgsq(N,A,b,u,EPS);
  else if ( BICGSTAB )
    its = bicgstab(N,A,b,u,EPS);
  else if ( GMRES )
    its = gmres(M,N,A,b,u,EPS);        // f"ur GMRES keine Vorkonditionierung
  else {
    cerr<<"kein L\"oser (-cghs -bicgsq -bicgstab -gmres) angegeben"<<endl;
    exit(-1);
  }

  {  // Residuum
    mult(A,u,v);
    for ( int i=0; i<N; ++i )
      w[i] = v[i]-b[i];
  }

  double c=0., d=0.;
  { // L2-Norm des Residuums
    for ( int i=0; i<N; ++i )
      c += w[i]*w[i];
    for ( int j=0; j<N; ++j )
      d += b[j]*b[j];
  }
  cout<<its<<" Iterationen, Residuum="<<sqrt(c)<<" Norm(b)="<<sqrt(d)<<endl;

  float used = (float)(clock() - start)/(float)(CLOCKS_PER_SEC);
  cout<<"ben\"otigte Zeit: "<<used<<endl;
  delete[] u;
  delete[] v;
  delete[] w;
  delete[] b;
  destructAntiMat(A);
}
