
NRGANDDIR=$(HOME)/Cprogs/NRG/Anderson/
CLIBDIR=$(HOME)/Cprogs/Clibrary/
ANDLIBDIR=$(NRGANDDIR)lib/
DATSDIR=$(NRGANDDIR)Dats/

TWOCHDIR=$(HOME)/Cprogs/TwoChNRG/
LAPACKDIR=/usr/src/lapack-3.0-5/lapack-3.0/SRC/


##OBJS_PREV=ProcReadParms.o \
##     ReadMatVec.o ReadQSParameters.o ReadQSVectors.o \
##     AllocThings.o ReallocThings.o


OBJS=NRGarrayClass.o NRGmatrixClass.o ij2r.o r2ij.o \
     CutStates.o QS_BuildBasis.o \
     dEqual.o CalcThermo.o \
     Sdots.o QfdQm1.o BuildBasis.o CGordan.o \
     UpdateMatrices.o FindMatchBlock.o OpMatRules.o \
     Calc_phMatEl.o DiagHN.o OpTables.o \
     FileExists.o NRGthermoClass.o


OBJS_TwoChQS=TwoChQS_SetSingleSite.o TwoChQS_SetH0.o \
     TwoChQS_UpdateMatrices.o TwoChQS_DiagHN.o \
     TwoChQS_OpMatRules.o TwoChQS_OpTables.o

CXXFLAGS = -g -I$(CLIBDIR) -I$(TWOCHDIR)src/ -I$(TWOCHDIR)include/ -I$(LAPACKDIR)
LDFLAGS = -lgsl -lgslcblas -llapack -lblas -lm
##CXX=gcc4
##LDFLAGS = -lstdc++ -lgsl -lgslcblas -llapack -lblas -lg2c -lm
## -loption SHOULD be added AT THE END!!! Stupid Make... 

##
## Useful: -v option in g++
##
## SOME IMPROVEMENTS to make things easier:
##
#
## Specify Directories that make will search

##VPATH=$(CLIBDIR):$(ANDLIBDIR):$(DATSDIR)
VPATH=$(CLIBDIR):$(TWOCHDIR)src/:$(TWOCHDIR)include/

##
## New pattern rule for "c" files
##
%.o : %.cpp
##	$(CXX) -x c++ -c -I./ $(CXXFLAGS) $(LDFLAGS) $< -o $@
	$(CXX) -c -I./ $(CXXFLAGS) $< -o $@
##
## $^  all dependencies
## $<  just the first dependency
## $@  name of the target



##OneChQS : OneChQS.cpp $(OBJS) $(OBJS_OneChQS)
##	$(CXX)  $^ $(CXXFLAGS) $(LDFLAGS) -o $@ 
##
## Don't know why it doesn't put LDFLAGS after the objects
##

TwoChQS : TwoChQS.cpp $(OBJS) $(OBJS_TwoChQS)
	$(CXX)  $^ $(CXXFLAGS) $(LDFLAGS) -o $@ 

TwoChQS_test : TwoChQS.cpp $(OBJS) $(OBJS_TwoChQS)
	$(CXX)  $^ $(CXXFLAGS) $(LDFLAGS) -o $@ 

clean :
	echo $(DATSDIR)
	-rm *.o *~
