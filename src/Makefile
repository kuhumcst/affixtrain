VERSION=1
MINOR_NUMBER=1

BIN=../

AFFIXTRAIN=$(BIN)affixtrain
AFFIXTRAINPROFILING=$(BIN)affixtrainp

CC=g++ -O3 -Wall -DNDEBUG -pedantic -std=c++11
CCPROFILING=g++ -O3 -Wall -c -pg -DNDEBUG -pedantic -std=c++11

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
	hashtable.cpp\
	$(LETTERFUNCDIR)/isofunc.cpp\
	$(LETTERFUNCDIR)/letter.cpp\
	$(LETTERFUNCDIR)/letterfunc.cpp\
	node.cpp\
	optionaff.cpp\
	ruletemplate.cpp\
	shortrulepair.cpp\
	strng.cpp\
	testrules.cpp\
	trainingpair.cpp\
	$(LETTERFUNCDIR)/utf8func.cpp\
	vertex.cpp

MAKEAFFIXRULESHEADERS=\
	affixtrain.h\
	applyaffrules.h\
	argopt.h\
	comp.h\
	flexcombi.h\
	hashtable.h\
	$(LETTERFUNCDIR)/isofunc.h\
	$(LETTERFUNCDIR)/letter.h\
	$(LETTERFUNCDIR)/letterfunc.h\
	node.h\
	optionaff.h\
	ruletemplate.h\
	settingsaffixtrain.h\
	shortrulepair.h\
	strng.h\
	testrules.h\
	trainingpair.h\
	$(LETTERFUNCDIR)/utf8func.h\
	vertex.h

MAKEAFFIXRULESOBJ=\
	affixtrain.o\
	applyaffrules.o\
	argopt.o\
	comp.o\
	flexcombi.o\
	hashtable.o\
	isofunc.o\
	letter.o\
	letterfunc.o\
	node.o\
	optionaff.o\
	ruletemplate.o\
	shortrulepair.o\
	strng.o\
	testrules.o\
	trainingpair.o\
	utf8func.o\
	vertex.o
	
				

#                               #
#   Dependencies start here.    #
#                               #


$(AFFIXTRAIN): $(MAKEAFFIXRULESSRC) $(MAKEAFFIXRULESHEADERS)
	$(CC) $(PIC) $(DEBUG) $(GCCINC) -c $(MAKEAFFIXRULESSRC)
	$(CCLINKDYNAMIC) $(MAKEAFFIXRULESOBJ) -o $@ $(GCCLINK)

$(AFFIXTRAINPROFILING): $(MAKEAFFIXRULESSRC) $(MAKEAFFIXRULESHEADERS)
	$(CCPROFILING) $(PIC) $(DEBUG) $(GCCINC) -c $(MAKEAFFIXRULESSRC)
	$(CCLINKDYNAMICP) $(MAKEAFFIXRULESOBJ) -o $@ $(GCCLINK)

all: $(AFFIXTRAIN) $(AFFIXTRAINPROFILING)

clean:
	$(RM) *.o
	$(RM) $(AFFIXTRAIN) $(AFFIXTRAINPROFILING)
