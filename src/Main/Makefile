
#NRGANDDIR=$(HOME)/Cprogs/NRG/Anderson/
#CLIBDIR=$(HOME)/Cprogs/Clibrary/
#ANDLIBDIR=$(NRGANDDIR)lib/
#DATSDIR=$(NRGANDDIR)Dats/

TWOCHDIR=$(HOME)/Cprogs/TwoChNRG/
##LAPACKDIR=/usr/src/lapack-3.0-5/lapack-3.0/SRC/


##OBJS_PREV=ProcReadParms.o \
##     ReadMatVec.o ReadQSParameters.o ReadQSVectors.o \
##     AllocThings.o ReallocThings.o


OBJS=NRGarrayClass.o NRGmatrixClass.o NRGthermoClass.o NRGchainClass.o \
     ij2r.o r2ij.o \
     CutStates.o QS_BuildBasis.o CGordan_SStot.o CGordan.o \
     dEqual.o  CalcThermo.o Sdots.o QfdQm1.o OpTables.o \
     UpdateMatrices.o FindMatchBlock.o CalcOpAvg.o FileExists.o \
     DiagHN.o BuildBasis.o NRGcodehandlerClass.o \
     OpMatRules.o Calc_phMatEl.o \
     OneChQS_SetSingleSite.o \
     TwoChQS_SetSingleSite.o TwoChQS_OpTables.o TwoChQS_OpMatRules.o \
     CalcDN.o CodeHandler_ModelSwitch.o HybFuncDefs.o SpecFuncClass.o \
     GSL_Integrator.o SpecFunc_CFSmethod.o SpecFunc_allSTL.o RotateMatrix.o 

OBJS_MAIN=Main_SetH0.o Main_OpMatRules.o Main_SetSingleSite.o

#CXXFLAGS =  -I$(CLIBDIR) -I$(TWOCHDIR)src/ -I$(TWOCHDIR)include/ \
#            -I$(TWOCHDIR)src/OneChQS/ -I$(TWOCHDIR)src/TwoChQS/ \
#            -I$(LAPACKDIR) -funroll-loops -O2
CXXFLAGS =  -I$(TWOCHDIR)include/ -I$(TWOCHDIR)src/ \
            -funroll-loops -O2
##LDFLAGS =  -lgsl -lgslcblas -llapack -lblas -lg2c -lm
# Ubuntu 8.10 does NOT have g77 (or g2c)
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
##VPATH=$(CLIBDIR):$(TWOCHDIR)src/:$(TWOCHDIR)include/:$(TWOCHDIR)src/OneChQS/:$(TWOCHDIR)src/TwoChQS/
## Removing OneChQS  and TwoChQS dirs (cpp files should be there...)
VPATH=$(CLIBDIR):$(TWOCHDIR)src/:$(TWOCHDIR)include/

##
## New pattern rule for "c" files
##
#%.o : %.c
#	$(CXX) -c -I./  $(CXXFLAGS) $(LDFLAGS) $< -o $@
%.o : %.cpp
	$(CXX) -c   $(CXXFLAGS) $< -o $@
##
## $^  all dependencies
## $<  just the first dependency
## $@  name of the target



NRG_main : NRG_main.cpp $(OBJS) $(OBJS_MAIN)
	$(CXX)  $^ $(CXXFLAGS) $(LDFLAGS) -o $@ 
##
## Don't know why it doesn't put LDFLAGS after the objects
##


NRG_main_bkup : NRG_main_bkup.cpp $(OBJS) $(OBJS_MAIN)
	$(CXX)  $^ $(CXXFLAGS) $(LDFLAGS) -o $@ 

clean :
	-rm *.o *~
