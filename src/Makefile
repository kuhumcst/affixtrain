VERSION=1
MINOR_NUMBER=1

BIN=~/bin/
#BIN=

AFFIXTRAIN=$(BIN)affixtrain

CC=g++ -O3 -Wall -DNDEBUG -pedantic

CCLINKDYNAMIC=g++ 

#DEBUG=-g
DEBUG=

GCCLINK=-lstdc++

LETTERFUNCDIR = ../../letterfunc/src
GCCINC += -I$(LETTERFUNCDIR)

RM=rm -f

				
MAKEAFFIXRULESSRC=\
	applyaffrules.cpp\
	argopt.cpp\
	comp.cpp\
	flexcombi.cpp\
	graph.cpp\
	$(LETTERFUNCDIR)/isofunc.cpp\
	$(LETTERFUNCDIR)/letter.cpp\
	$(LETTERFUNCDIR)/letterfunc.cpp\
	optionaff.cpp\
	affixtrain.cpp\
	$(LETTERFUNCDIR)/utf8func.cpp

MAKEAFFIXRULESOBJ=\
	applyaffrules.o\
	argopt.o\
	comp.o\
	graph.o\
	flexcombi.o\
	isofunc.o\
	letter.o\
	letterfunc.o\
	optionaff.o\
	affixtrain.o\
	utf8func.o
				

#                               #
#   Dependencies start here.    #
#                               #


$(AFFIXTRAIN): $(MAKEAFFIXRULESSRC)
	$(CC) $(PIC) $(DEBUG) $(GCCINC) -c $(MAKEAFFIXRULESSRC)
	$(CCLINKDYNAMIC) $(MAKEAFFIXRULESOBJ) -o $@ $(GCCLINK)

all: $(AFFIXTRAIN)

clean:
	$(RM) *.o
	$(RM) $(AFFIXTRAIN)
