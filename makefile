#
# This is the makefile for libcash2-s.so
#

###############################################
# Specify where you want to install the library
#INSTALL_DIRECTORY = ../cash

##############################
# Specify which program to use
CPP = g++
PPOPT = -w -ggdb3 -O3  #Wall

CC = gcc
CCOPT = -c -O2 #-Wall #-fPIC

SED = sed
INSTALL = install
RM = /bin/rm -f
RMDIR = rmdir
#########################################################
# Do not modify below unless you know what you want to do
OBJCASH = $(BREM_CASH)arithmetic.o $(BREM_CASH)basic.o $(BREM_CASH)color.o $(BREM_CASH)filter.o $(BREM_CASH)io.o $(BREM_CASH)logical.o \
$(BREM_CASH)margolus.o $(BREM_CASH)movie.o $(BREM_CASH)neighbors.o $(BREM_CASH)noise.o $(BREM_CASH)png.o $(BREM_CASH)ps.o $(BREM_CASH)random.o $(BREM_CASH)shift.o\
$(BREM_CASH)x11.o
OBJPOAS = $(POAS)/Genome.o $(POAS)/Pearl.o $(POAS)/HK.o $(POAS)/NonEss.o $(POAS)/eDNA.o $(POAS)/Noncoding.o $(POAS)/Transposase.o
OBJCASH2 = $(BREM_CASH)cash2.o $(BREM_CASH)mersenne.o
OBJMAIN = $(BREM_CASH)cash2-s.o

MYOBJ=selfish.cpp clibs.hpp functions.hpp functions_data.hpp gotmouse.cpp

LDIR = -L/usr/X11R6/lib
LIBS = -lpng -lz -lX11 -lm -lgrace_np -lm

BREM_CASH = Brem-cash/
POAS = PoaS/

OBJVIB = selfish.o

all: selfish

rule: dependencies
	command
	other command
	@silent command

selfish.o:  ${MYOBJ}
	$(CPP) $(PPOPT) -c selfish.cpp

$(BREM_CASH)%.o: $(BREM_CASH)%.c
	$(CC) $(CCOPT) -o $@ $<

$(BREM_CASH)%.o: $(BREM_CASH)%.cpp $(BREM_CASH)%.cc
	$(CPP) -c $(PPOPT) -o $@ $<

$(POAS)%.o: $(POAS)%.cc
	$(CPP) -c $(PPOPT) -o $@ $<

$(OBJCASH): $(BREM_CASH)cash2003.h
$(OBJCASH2): $(BREM_CASH)cash2.hpp $(BREM_CASH)mersenne.h
$(OBJMAIN): $(BREM_CASH)cash2003.h $(BREM_CASH)cash2.hpp $(BREM_CASH)mersenne.h $(BREM_CASH)cash2-s.hpp


selfish: $(OBJCASH) $(OBJCASH2) $(OBJMAIN) $(OBJVIB) $(OBJPOAS)
	$(CPP) $(PPOPT)  $(OBJVIB) $(OBJCASH) $(OBJCASH2) $(OBJMAIN) $(OBJPOAS) -o SelfishDNA \
	  -lm  -lc -lpng -lz -lX11 -lgrace_np -lcurses

.PHONY: selfish.cpp
