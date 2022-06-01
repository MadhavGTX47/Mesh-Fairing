//
// 3dim Laplace mit lin Elementen
//

#include <assert.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
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



// ============================================================================

template< class t > inline t max( t a, t b ) { return (a>b) ? a : b; }

// ============================================================================


const double d23=2./3., d16=-1./6., d00=0.;
const double d827=8./27., d427=4./27., d227=2./27., d127=1./27.;

const double SM[8][8] = { {d23,d00,d16,d00,d00,d16,d16,d16},
			  {d00,d23,d00,d16,d16,d00,d16,d16},
			  {d16,d00,d23,d00,d16,d16,d00,d16},
			  {d00,d16,d00,d23,d16,d16,d16,d00},
			  {d00,d16,d16,d16,d23,d00,d16,d00},
			  {d16,d00,d16,d16,d00,d23,d00,d16},
			  {d16,d16,d00,d16,d16,d00,d23,d00},
			  {d16,d16,d16,d00,d00,d16,d00,d23} };

const double MM[8][8] = { {d827,d427,d227,d427,d427,d227,d127,d227},
			  {d427,d827,d427,d227,d227,d427,d227,d127},
			  {d227,d427,d827,d427,d127,d227,d427,d227},
			  {d427,d227,d427,d827,d227,d127,d227,d427},
			  {d427,d227,d127,d227,d827,d427,d227,d427},
			  {d227,d427,d227,d127,d427,d827,d427,d227},
			  {d127,d227,d427,d227,d227,d427,d827,d427},
			  {d227,d127,d227,d427,d427,d227,d427,d827} };

void mult( const double (&A)[8][8], double *v, double *w ) {
  for ( int i=0; i<8; ++i ) {
    double tmp=0.;
    for ( int j=0; j<8; ++j )
      tmp+=A[i][j]*v[j];
    w[i]=tmp;
  }
}

const int ISF[8][3]    = { {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
			   {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1} };

// ============================================================================


struct Element {
  int v[8];                     // Index der 8 Vertices
  double reg1[8], reg2[8];      // Zwischenspeicher f"ur Daten
};



// ============================================================================


struct Matrix0bnd {
  Matrix0bnd( const double (&M)[8][8], double h,
	      int NV, int NE, Element *E, int NB, int *B )
    : M(M), NE(NE), NV(NV), NB(NB), E(E), B(B), h(h) { }
  const double (&M)[8][8];      // Element-Matrix
  int NE, NV, NB;               // #Elementen, #Vertices, #Randvertices
  Element *E;                   // Elemente
  int *B;                       // Indices der Randvertices
  double h;                     // Verkleinerungs-Faktor
};

void mult( const Matrix0bnd &A, double *v, double *w ) {
  { for ( int k=0; k<A.NB; ++k ) v[A.B[k]]=0.; }        // Nullrand v
  { for ( int k=0; k<A.NV; ++k ) w[k]=0.; }             // w=0
  { for ( int k=0; k<A.NE; ++k ) {
    Element &E(A.E[k]);
    { for ( int i=0; i<8; ++i )  E.reg1[i]=v[E.v[i]]; } // Vektor=>Element
    mult(A.M,E.reg1,E.reg2);                            // Elementmatrix
    { for ( int i=0; i<8; ++i )  w[E.v[i]]+=E.reg2[i]; }// Element=>Vektor
  } }
  { for ( int k=0; k<A.NB; ++k ) w[A.B[k]]=0.; }        // Nullrand w
  { for ( int k=0; k<A.NV; ++k ) w[k]*=A.h;    }        // Faktor h
}


// ----------------------------------------------------------------------------

struct StiffMatrix
  : Matrix0bnd {
  StiffMatrix( double h, int NV, int NE, Element *E, int NB, int *B )
    : Matrix0bnd(SM,h,NV,NE,E,NB,B) { };
};

struct MassMatrix
  : Matrix0bnd {
  MassMatrix( double h, int NV, int NE, Element *E, int NB, int *B )
    : Matrix0bnd(MM,h*h*h,NV,NE,E,NB,B) { };
};

// ============================================================================


inline int posV( int n, int i, int j, int k ) {
  return (n+1)*( (n+1)*i + j ) + k;
}

inline int posE( int n, int i, int j, int k ) {
  return n*( n*i + j ) + k;
}


// ============================================================================

typedef double double3[3];
inline void setdouble3( double3 &D, double d1, double d2, double d3 ) {
  D[0]=d1, D[1]=d2, D[2]=d3;
}

// ============================================================================


const double pi = 3.14159265358979323846;

#if 0
// auf [-1,1]^3
struct SinSinSinSolution {
  double operator()( double x, double y, double z ) const {
    return sin(pi*x)*sin(pi*y)*sin(pi*z);
  }
} g;
struct SinSinSin {
  double operator()( double x, double y, double z ) const {
    return 3.*pi*pi*sin(pi*x)*sin(pi*y)*sin(pi*z);
  }
} f;
#else
// auf [-1,1]^3
struct CosCosCosSolution {
  double operator()( double x, double y, double z ) const {
    return cos(.5*pi*x)*cos(.5*pi*y)*cos(.5*pi*z);
  }
} g;
struct CosCosCos {
  double operator()( double x, double y, double z ) const {
    return .25*3.*pi*pi*cos(.5*pi*x)*cos(.5*pi*y)*cos(.5*pi*z);
  }
} f;
#endif


int main( int argc, char *argv[] ) {
  Arguments args(argc,argv);
  int n         = args.get_int_argument("-n",50);  // #Punkte in eine Richtung
  int m         = args.get_int_argument("-m",5);   // f"ur gmres(M)
  bool DETAILED = args.contains_flag("-detailed"); // ausf"uhrliche Iterationen
  bool CGHS     = args.contains_flag("-cghs");
  bool BICGSQ   = args.contains_flag("-bicgsq");
  bool BICGSTAB = args.contains_flag("-bicgstab");
  bool GMRES    = args.contains_flag("-gmres");
  double EPS    = args.get_double_argument("-eps",1e-5);
  if ( !args.no_more_args() ) {
    cerr<<"falsche Argumente"<<endl;
    exit(-1);
  }

  clock_t t=clock();              // Anfangszeitpunkt
  int NV=(n+1)*(n+1)*(n+1);       // #Vertices
  int NE=n*n*n;                   // #Elemente
  int NB=6*(n+1)*(n+1);           // #Randvertices
  double h=1./double(n);          // Gitterweite=2*h

  double3 *V=0;                   // Koord der Vertices
  int     *B=0;                   // Indices der Randpunkte
  Element *E=0;                   // die Elemente

  double *c=0;                    // exakte L"osung
  double *b=0;                    // rechte Seite
  double *x=0;                    // L"osungsvektor

  int its=-1;                     // Iterationen L"oser
  double diff=0.;                 // Abweichung von exakter L"osung
  double time=0.;                 // ben"otigte Zeit in Sekunden

  { // Vertex-Koordinaten setzen
    V = new double3[NV];
    for ( int i=0; i<=n; ++i )
      for ( int j=0; j<=n; ++j )
	for ( int k=0; k<=n; ++k )
	  setdouble3( V[posV(n,i,j,k)],
		      -1.+2.*double(i)*h,
		      -1.+2.*double(j)*h,
		      -1.+2.*double(k)*h );
  }
  { // Elemente setzen
    E = new Element[NE];
    for ( int i=0; i<n; ++i )
      for ( int j=0; j<n; ++j )
	for ( int k=0; k<n; ++k )
	  for ( int v=0; v<8; ++v )
	    E[posE(n,i,j,k)].v[v]
	      = posV(n,i+ISF[v][0],j+ISF[v][1],k+ISF[v][2]);
  }
  { // Rand setzen (die Punkte an den Kanten und Ecken mehrfach - macht nichts)
    B = new int[NB];
    int k=0;
    for ( int i=0; i<=n; ++i )
      for ( int j=0; j<=n; ++j )
	{
	  B[k++] = posV(n,0,i,j);
	  B[k++] = posV(n,n,i,j);
	  B[k++] = posV(n,i,0,j);
	  B[k++] = posV(n,i,n,j);
	  B[k++] = posV(n,i,j,0);
	  B[k++] = posV(n,i,j,n);
	}
  }

  { // rechte Seite setzen
    double *aux = new double[NV];
    for ( int i=0; i<=n; ++i )
      for ( int j=0; j<=n; ++j )
	for ( int k=0; k<=n; ++k )
	  {
	    int m = posV(n,i,j,k);
	    aux[m] = f(V[m][0],V[m][1],V[m][2]);
	  }
    MassMatrix MM(h,NV,NE,E,NB,B);
    b = new double[NV];
    mult(MM,aux,b);                   // Integrale = ca. Mult mit Massenmatrix
  }

  {
    // L"osungsvektor initialisieren
    x = new double[NV];
    for ( int k=0; k<NV; ++k )
      x[k] = 0.;
  }

  { // und l"osen
    StiffMatrix SM(h,NV,NE,E,NB,B);
    if ( CGHS )
      its = cghs(NV,SM,b,x,EPS,DETAILED);
    else if ( BICGSQ )
      its = bicgsq(NV,SM,b,x,EPS,DETAILED);
    else if ( BICGSTAB )
      its = bicgstab(NV,SM,b,x,EPS,DETAILED);
    else if ( GMRES )
      its = gmres(m,NV,SM,b,x,EPS,DETAILED);  // GMRES ohne Vorkonditionierer
    else {
      cerr<<"kein L\"oser (-cghs -bicgsq -bicgstab -gmres) angegeben"<<endl;
      exit(-1);
    }
  }

  { // exakte L"osung setzen
    c = new double[NV];
    for ( int i=0; i<=n; ++i )
      for ( int j=0; j<=n; ++j )
	for ( int k=0; k<=n; ++k )
	  {
	    int m = posV(n,i,j,k);
	    c[m] = g(V[m][0],V[m][1],V[m][2]);
	  }
  }

#if 0
  { // alle Daten ausgeben
    cout<<setprecision(5);
    for ( int k=0; k<NV; ++k )
      cout<<setw(20)<<b[k]
	  <<setw(20)<<c[k]
	  <<setw(20)<<x[k]
	  <<endl;
  }
#endif

  { // und Abweichung bestimmen
    diff=0.;
    for ( int k=0; k<NV; ++k )
      diff = max( diff, fabs(c[k]-x[k]) );
  }

  { // ben"otigte Zeit berechnen
    time = double(clock()-t)/double(CLOCKS_PER_SEC);
  }

  { // und alles ausgeben
    cout<<"n="<<n
	<<" \tdofs="<<NV
	<<" \tits="<<its
	<<" \tdiff="<<diff
	<<" \th="<<2.*h
	<<" \ttime="<<time
	<<endl;
  }
  { // wieder alles deleten
    delete[] V;
    delete[] E;
    delete[] B;
    delete[] b;
    delete[] x;
    delete[] c;
  }
}
