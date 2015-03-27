VERSION=1
MINOR_NUMBER=1

BIN=~/bin/
#BIN=

AFFIXTRAIN=$(BIN)affixtrain
AFFIXTRAINI=$(BIN)affixtraini
AFFIXTRAINP=$(BIN)affixtrainp

CC=g++ -O3 -Wall -DNDEBUG -pedantic
CCP=g++ -O3 -Wall -c -pg -DNDEBUG -pedantic

CCLINKDYNAMIC=g++ 
CCLINKDYNAMICP=g++ -Wall -pg 

#DEBUG=-g
DEBUG=

GCCLINK=-lstdc++

LETTERFUNCDIR = ../../letterfunc/src
GCCINC += -I$(LETTERFUNCDIR)

RM=rm -f

				
MAKEAFFIXRULESSRC=\
	affixtrain.cpp\
	applyaffrules.cpp\
	argopt.cpp\
	comp.cpp\
	flexcombi.cpp\
	graph.cpp\
	testrules.cpp\
	$(LETTERFUNCDIR)/isofunc.cpp\
	$(LETTERFUNCDIR)/letter.cpp\
	$(LETTERFUNCDIR)/letterfunc.cpp\
	optionaff.cpp\
	$(LETTERFUNCDIR)/utf8func.cpp

MAKEAFFIXRULESHEADERS=\
	applyaffrules.h\
	argopt.h\
	comp.h\
	flexcombi.h\
	graph.h\
	affixtrain.h\
	testrules.h\
	$(LETTERFUNCDIR)/isofunc.h\
	$(LETTERFUNCDIR)/letter.h\
	$(LETTERFUNCDIR)/letterfunc.h\
	optionaff.h\
	settingsaffixtrain.h\
	strng.h\
	$(LETTERFUNCDIR)/utf8func.h

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
	testrules.o\
	utf8func.o
				

#                               #
#   Dependencies start here.    #
#                               #


$(AFFIXTRAIN): $(MAKEAFFIXRULESSRC) $(MAKEAFFIXRULESHEADERS)
	$(CC) $(PIC) $(DEBUG) $(GCCINC) -c $(MAKEAFFIXRULESSRC)
	$(CCLINKDYNAMIC) $(MAKEAFFIXRULESOBJ) -o $@ $(GCCLINK)

$(AFFIXTRAINI): $(MAKEAFFIXRULESSRC) $(MAKEAFFIXRULESHEADERS)
	$(CC) $(PIC) $(DEBUG) $(GCCINC) -c $(MAKEAFFIXRULESSRC) -D INTPARMS
	$(CCLINKDYNAMIC) $(MAKEAFFIXRULESOBJ) -o $@ $(GCCLINK)

$(AFFIXTRAINP): $(MAKEAFFIXRULESSRC) $(MAKEAFFIXRULESHEADERS)
	$(CCP) $(PIC) $(DEBUG) $(GCCINC) -c $(MAKEAFFIXRULESSRC)
	$(CCLINKDYNAMICP) $(MAKEAFFIXRULESOBJ) -o $@ $(GCCLINK)

all: $(AFFIXTRAIN) $(AFFIXTRAINP)

clean:
	$(RM) *.o
	$(RM) $(AFFIXTRAIN) $(AFFIXTRAINI)
