
NRGANDDIR=$(HOME)/Cprogs/NRG/Anderson/
CLIBDIR=$(HOME)/Cprogs/Clibrary/
ANDLIBDIR=$(NRGANDDIR)lib/
DATSDIR=$(NRGANDDIR)Dats/

TWOCHDIR=$(HOME)/Cprogs/TwoChNRG/
LAPACKDIR=/usr/src/lapack-3.0-5/lapack-3.0/SRC/


##OBJS_PREV=ProcReadParms.o \
##     ReadMatVec.o ReadQSParameters.o ReadQSVectors.o \
##     AllocThings.o ReallocThings.o


OBJS=NRGarrayClass.o NRGmatrixClass.o NRGthermoClass.o NRGchainClass.o \
     ij2r.o r2ij.o \
     CutStates.o QS_BuildBasis.o CGordan_SStot.o CGordan.o \
     dEqual.o  CalcThermo.o Sdots.o QfdQm1.o OpTables.o \
     UpdateMatrices.o FindMatchBlock.o CalcOpAvg.o FileExists.o \
     DiagHN.o

OBJS_OneChQS=OneChQS_BuildBasis.o OneChQS_DiagHN.o \
     OneChQS_prefactorHN.o OneChQS_UpdateMatrices.o \
     OneChQS_SetSingleSite.o OneChQS_OpTables.o \
     OneChQS_SetH0.o OneChQS_OpMatRules.o

CXXFLAGS =  -I$(CLIBDIR) -I$(TWOCHDIR)src/ -I$(TWOCHDIR)include/ \
            -I$(LAPACKDIR) -funroll-loops -O2
LDFLAGS =  -lgsl -lgslcblas -llapack -lblas -lm
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
%.o : %.c
	$(CXX) -c -I./ $(CXXFLAGS) $(LDFLAGS) $< -o $@
##
## $^  all dependencies
## $<  just the first dependency
## $@  name of the target



OneChQS : OneChQS.cpp $(OBJS) $(OBJS_OneChQS)
	$(CXX)  $^ $(CXXFLAGS) $(LDFLAGS) -o $@ 
##
## Don't know why it doesn't put LDFLAGS after the objects
##

Kondo_OneChQS : Kondo_OneChQS.cpp $(OBJS) $(OBJS_OneChQS)
	$(CXX)  $^ $(CXXFLAGS) $(LDFLAGS) -o $@ 


clean :
	echo $(DATSDIR)
	-rm *.o *~
