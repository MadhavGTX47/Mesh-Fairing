###############################################################################
#
#                dieses Makefile ben"otigt GNU make!!!
#
###############################################################################
#
# wenn auf der SGI kompiliert wird: Wahl zwischen GNU g++ und SGI-CC
#USE_COMPILER = GNU
USE_COMPILER = SGI
#
# optimiert kompilieren oder zum Debuggen
#FLAGTYPE = OPT
FLAGTYPE = DEBUG
#
###############################################################################

  CXX      = g++ -W -ansi -pedantic -Wall 
  OPTFLAGS = -O3 -funroll-loops -fast
  DBGFLAGS = -g
  FLIB     = -framework vecLib
  BLAS	   = -lblas
  DEPFLAG  = -MM




ifeq ($(FLAGTYPE),OPT)
  CXXFLAGS = $(OPTFLAGS)
endif

ifeq ($(FLAGTYPE),DEBUG)
  CXXFLAGS = $(DBGFLAGS)
endif

MKDEP = $(CXX) $(DEPFLAG)

###############################################################################



# hier werden die Dependencies gespeichert
DEPFILE	 = Dependfile

# Name, der beim tarren verwendet wird (ls = linear solver)
PCKNAME  = ls

# heutiges Datum
DATE     = $(shell date +%D | sed 's/\// /g' | awk '{print $$2$$1$$3}')

# Files, die getarred werden sollen
TARFILES = README *.h www example

# Name des Tarfiles
TARNAME	 = $(PCKNAME)$(DATE)

# alle Sourcen
SRC	 = antisym.cc laplace2d.cc laplace3d.cc argument.cc

# make all
all: depend antisym laplace2d laplace3d


laplace3d: laplace3d.o argument.o
	$(CXX) $(CXXFLAGS) -o laplace3d laplace3d.o argument.o \
	-lm $(BLAS) $(FLIB)

laplace2d: laplace2d.o argument.o
	$(CXX) $(CXXFLAGS) -o laplace2d laplace2d.o argument.o \
	-lm $(BLAS) $(FLIB)

antisym: antisym.o argument.o
	$(CXX) $(CXXFLAGS) -o antisym antisym.o argument.o \
	-lm $(BLAS) $(FLIB)

# tgz-File erzeugen (tarren und gzippen)
tgz:
	cd ..; tar -cf - $(TARFILES) | gzip >$(TARNAME).tgz;

# tgz-File erzeugen und symbolic link vom www-Verzeichnis darauf setzen
dist: tgz
	cd ../www; rm -f lsolver.tgz; ln -s ../$(TARNAME).tgz lsolver.tgz

# Dependencies erstellen
depend:
	$(MKDEP) $(SRC) > $(DEPFILE)

# aufr"aumen
clean:
	rm -f *.o *~ core

# mehr aufr"aumen
proper: clean
	rm -rf ii_files $(DEPFILE) antisym laplace2d laplace3d

# Dependencies einf"ugen
-include $(DEPFILE)
