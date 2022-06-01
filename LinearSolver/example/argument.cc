#include <cstdlib>
#include <string>
#include <assert.h>
#include <vector>
#include "argument.h"

#ifndef __KCC
#  include <algo.h>
#else
#  include <algobase.h>
#  include <algorithm.h>
#  include <functional.h>
#endif


// ============================================================================


struct ArgumentsData {
  typedef const char *string;
  ArgumentsData( int argc, char *argv[] );
  ArgumentsData( void );
  int             argc;
  const string   *cargv;
	std::vector<string>  argv;
  std::vector<string*> aux;
};

struct string_eq {
  typedef const char *string;
  string_eq( string s ) : s(s) { }
  bool operator()( string t ) const { return t!=0 && 0==strcmp(s,t); }
  string s;
};


// ============================================================================


ArgumentsData::
ArgumentsData( int argc, char *argv[] )
  : argc(argc-1), cargv(argv), argv(argc-1), aux(0)
{
  //  copy(cargv+1,cargv+1+Arguments::argc,Arguments::argv);
  for ( int i=0; i<ArgumentsData::argc; ++i )
    ArgumentsData::argv[i]=argv[i+1];
}


ArgumentsData::
ArgumentsData( void ) {
  for ( std::vector<string*>::iterator pos=aux.begin(); pos!=aux.end(); ++pos )
    delete[] *pos;
}


// ============================================================================


Arguments::
Arguments( int argc, char *argv[] ) {
  data=new ArgumentsData(argc,argv);
}


Arguments::
~Arguments( void ) {
  delete data;
  data=0;
}


const char *Arguments::
get_program_name( void ) const {
  return data->cargv[0];
}


long Arguments::
get_long_argument( const char *argname, long default_value ) const {
	std::vector<string> &argv((std::vector<string>&)data->argv);
	std::vector<string>::iterator
    npos=find_if( argv.begin(), argv.end(), string_eq(argname) ),
    apos=npos;
  ++apos;
  char *ptr=0;
  long ret=default_value;
  if ( npos!=argv.end() && apos!=argv.end() ) {
    ret=strtol(*apos,&ptr,0);
    if ( ptr!=*apos )
      *npos=*apos=0;
  }
  return ret;
}


int Arguments::
get_int_argument( const char *argname, int default_value ) const {
  return (int)get_long_argument(argname,(long)default_value);
}


double Arguments::
get_double_argument( const char *argname, double default_value ) const {
	std::vector<string> &argv((std::vector<string>&)data->argv);
	std::vector<string>::iterator
    npos=find_if( argv.begin(), argv.end(), string_eq(argname) ),
    apos=npos;
  ++apos;
  char *ptr=0;
  double ret=default_value;
  if ( npos!=argv.end() && apos!=argv.end() ) {
    ret=strtod(*apos,&ptr);
    if ( ptr!=*apos )
      *npos=*apos=0;
  }
  return ret;
}


const char *Arguments::
get_string_argument( const char *argname, const char *default_value ) const {
	std::vector<string> &argv((std::vector<string>&)data->argv);
	std::vector<string>::iterator
    npos=find_if( argv.begin(), argv.end(), string_eq(argname) ),
    apos=npos;
  ++apos;
  string ret=default_value;
  if ( npos!=argv.end() && apos!=argv.end() )
    ret=*apos, *npos=*apos=0;
  return ret;
}


const char **Arguments::
get_n_string_arguments( int n, const char *argname ) const {
	std::vector<string> &argv((std::vector<string>&)data->argv);
	std::vector<string>::iterator
    npos=find_if( argv.begin(), argv.end(), string_eq(argname) ),
    apos=npos;
  int i=n;
  while ( i-->0 && ++npos!=argv.end() && *npos!=0 )
    ;
  string *ret;
  if ( i==0 )
    ret=(string*)data->cargv+2+(apos-argv.begin());
  else
    ret=0;
  return ret;
}


bool Arguments::
contains_flag( const char *flagname ) const {
	std::vector<string> &argv((std::vector<string>&)data->argv);
	std::vector<string>::iterator
    fpos=find_if( argv.begin(), argv.end(), string_eq(flagname) );
  bool ret=false;
  if ( fpos!=argv.end() && *fpos!=0 )
    ret=true, *fpos=0;
  return ret;
}


void Arguments::
get_rest( int &n, string *&rest ) const {
  const std::vector<string> &argv(data->argv);
	std::vector<string*>      &aux((std::vector<string*>&)data->aux);
#ifndef __KCC
  n=0;
  count_if( argv.begin(), argv.end(), bind2nd(not_equal_to<string>(),0), n );
#else
  n=count_if( argv.begin(), argv.end(), bind2nd(not_equal_to<string>(),0) );
#endif
  rest=new string[n];
  ((std::vector<string*>&)aux).push_back(rest);          // kein mutable vorhanden
  int k=0;
  for ( std::vector<string>::const_iterator pos=argv.begin();
	pos!=argv.end();
	++pos )
    if ( *pos!=0 )
      rest[k++]=*pos;
  assert(k==n);
}


bool Arguments::
no_more_args( void ) const {
  const std::vector<string> &argv(data->argv);
	std::vector<string>::const_iterator first_nonzero = 
    find_if( argv.begin(), argv.end(), bind2nd(not_equal_to<string>(),0) );
  return ( first_nonzero==argv.end() );
}


// ============================================================================
  
