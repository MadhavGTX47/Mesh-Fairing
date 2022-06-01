// Emacs soll diesen Header als -*- C++ -*- erkennen


#ifndef ARGUMENT_H
#define ARGUMENT_H
// ============================================================================
//
//  Arguments
//  =========
//
//  Einfacher Zugriff auf command line Argumente; falls ein Argument schon
//  abgefragt wurde, wird es aus der Liste gel"oscht, es steht dann nicht
//  mehr zur Verf"ugung.
//
//
//  Konstruktor:  Arguments( int argc, char *argv[] )
//     <argc> und <argv> sind die C-Argumente von   main(argc,argv)
//
//  string get_program_name( void ) const
//     Programmname (argv[0])
//
//  int get_int_argument( string argname, int default_value ) const
//     Wandle das nach <argname> folgende Argument in integer um; falls nicht
//     m"oglich (z.B. Argument existiert nicht), dann gebe <default_value>
//     zur"uck
//
//  long get_long_argument( string argname, long default_value ) const
//     entsprechend get_int_argument, nur long integer
//
//  double get_double_argument( string argname, double default_value ) const
//     entsprechend get_int_argument, jedoch double
//
//  string get_string_argument( string argname, string default_value ) const
//     entsprechend get_int_argument, jedoch string
//
//  string *get_n_string_arguments( int n, string argname ) const
//     Falls nach <argname> noch mind. <n> Argumente folgen, werden diese in
//     einem string-Array zur"uckgegeben (Teil-String von argv, bleibt also
//     bis zum Ende des Programms erhalten, darf aber nicht ge"andert werden).
//     Falls keine n Argumente mehr folgen, oder "uberhaupt kein Argument
//     <argname> auftritt, wird ein 0-Zeiger zur"uckgegeben.
//
//  bool contains_flag( string flagname ) const
//     true, falls das Flag <flagname> auftritt, false sonst
//
//  void get_rest( int &n, string *&rest ) const
//     <n> ent"halt die Anzahl der noch nicht abgefragten Argumente,
//     <rest> enth"alt String-Array auf diese Argumente (dieses Array bleibt
//     bis zum Zerst"oren von Arguments erhalten, wird dann jedoch
//     mitzerst"ort).
//
//  bool no_more_args( void ) const
//     true, falls schon alle Argumente abgefragt wurden
//
//
//  Christian Badura, 17.06.98
//
// ============================================================================


struct ArgumentsData;


class Arguments {

  typedef const char *string;

public:
  Arguments( int argc, char *argv[] );
  ~Arguments( void );

  string get_program_name( void ) const;
  int get_int_argument( string argname, int default_value ) const;
  long get_long_argument( string argname, long default_value ) const;
  double get_double_argument( string argname, double default_value ) const;
  string get_string_argument( string argname, string default_value ) const;
  string *get_n_string_arguments( int n, string argname ) const;
  bool contains_flag( string flagname ) const;
  void get_rest( int &n, string *&rest ) const;
  bool no_more_args( void ) const;

private:
  ArgumentsData  *data;

};


// ============================================================================
#endif // ARGUMENT_H
