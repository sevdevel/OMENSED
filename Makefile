# Generic Makefile for C-GEM
#
# Author - Sebastiaan van de Velde
# 
# Version - @{CGEM}Makefile   V0.1
#

# default machine dependent parameters.
FC =
FCDEFS =
FCDEBUG =
FCIFLAGS =
FLIB_FILES =
CPP =
CPPIFLAGS =
CPPF =
CPPOPTS =
CPPDEFS =

# options for cpp processing and libraries
include $(DOPTSFILE)
SOURCEMOD = main

# file paths
VPATH = source
CMPDIR = comps

# main program 
EXEFILE = main
MAINFILE = main_model.o

# no target
none:
	@echo No target provided

# invoke make
include $(CMPDIR)/compilers.cmp

# make command
MAKE = make

# object files
include $(CMPDIR)/objects.cmp
MODFILES = $(MODULES) 
SUBFILES = $(SUBPROGS) 

# compilation flags
CPPFLAGS = $(CPPIFLAGS) $(CPPOPTS) $(CPPDEFS)
FFLAGS = $(FCIFLAGS) $(FCDEFS) $(FCDEBUG)

# clean unnecessary files
clean:
	@rm -f $(MODFILES) $(SUBFILES) $(MAINFILE) $(EXEFILE) *.mod

# create executable
$(EXEFILE): $(MODFILES) $(SUBFILES) $(MAINFILE)
	$(FC) $(FFLAGS) -o $(EXEFILE) $(MODFILES) $(SUBFILES) \
	$(MAINFILE) $(FLIB_FILES)

# dependencies
include $(CMPDIR)/dependencies.cmp

# suffix rules
.SUFFIXES:
.SUFFIXES:.o .f90 .F90

.F90.f90:
	@rm -f $*.f90
	$(CPPF) $(CPPFLAGS) $< ./$*.f90
.f90.o:
	$(FC) $(FFLAGS) -c $< -o $*.o
.F90.o:
	@rm -f $*.f90
	$(CPPF) $(CPPFLAGS) $< ./$*.f90
	$(FC) $(FFLAGS) -c ./$*.f90 -o $*.o
	@rm -f $*.f90
