
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
     dEqual.o Q1Q2Sz_BuildBasis.o CalcSuscep.o \
     Sdots.o QfdQm1.o OpTables.o BuildBasis.o CGordan.o


OBJS_TwoChQ1Q2Sz=TwoChQ1Q2Sz_SetSingleSite.o 
	TwoChQ1Q2Sz_SetH0.o \
        TwoChQ1Q2Sz_OpTables.o 
	TwoChQ1Q2Sz_UpdateMatrices.o \
        TwoChQ1Q2Sz_DiagHN.o


CXXFLAGS = -I$(CLIBDIR) -I$(ANDLIBDIR) -I$(TWOCHDIR)include/ -I$(LAPACKDIR)
LDFLAGS =  -lgsl -lgslcblas -llapack -lblas -lg2c -lm
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



##OneChQS : OneChQS.cpp $(OBJS) $(OBJS_OneChQS)
##	$(CXX)  $^ $(CXXFLAGS) $(LDFLAGS) -o $@ 
##
## Don't know why it doesn't put LDFLAGS after the objects
##

TwoChQ1Q2Sz : TwoChQ1Q2Sz.cpp $(OBJS) $(OBJS_TwoChQ1Q2Sz)
	$(CXX)  $^ $(CXXFLAGS) $(LDFLAGS) -o $@ 


clean :
	echo $(DATSDIR)
	-rm *.o *~
