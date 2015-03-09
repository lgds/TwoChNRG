
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
     DiagHN.o BuildBasis.o NRGcodehandlerClass.o \
     OpMatRules.o  Calc_phMatEl.o \
     OneChQS_SetSingleSite.o \
     TwoChQS_SetSingleSite.o TwoChQS_OpTables.o TwoChQS_OpMatRules.o \
     SetSingleSite.o BroadDelta.o SpecFuncClass.o CalcDN.o RotateMatrix.o \
     HybFuncDefs.o GSL_Integrator.o SpecFunc_CFSmethod.o
##
##     TwoChQS_OpMatChecks.o # overlaps with TwoChQS_OpMatRules.o 

OBJS_DMNRG=DM_NRG_SetChildSt.o DM_NRG_CalcRhoN.o DM_NRG_SetRhoNmax.o \
            DM_NRG_CalcSpecFuncs.o DM_NRG_CommandLineRead.o 

OBJS_CONDUCTANCE=ConductanceClass.o

#CXXFLAGS =  -I$(CLIBDIR) -I$(TWOCHDIR)src/ -I$(TWOCHDIR)include/ \
#            -I$(TWOCHDIR)src/OneChQS/ -I$(TWOCHDIR)src/TwoChQS/ \
#            -I$(LAPACKDIR) -funroll-loops -O2
CXXFLAGS =  -I$(CLIBDIR) -I$(TWOCHDIR)include/ -I$(TWOCHDIR)src/ \
            -I$(LAPACKDIR) -funroll-loops -O2
## Do I really need -l g2c? Perhaps not (Ubuntu does not have it) 
##LDFLAGS =  -lgsl -lgslcblas -llapack -lblas -lg2c -lm
LDFLAGS =  -lgsl -lgslcblas -llapack -lblas  -lm
## -loption SHOULD be added AT THE END!!! Stupid Make... 


##
## Useful: -v option in g++
##
## SOME IMPROVEMENTS to make things easier:
##
#
## Specify Directories that make will search

##VPATH=$(CLIBDIR):$(ANDLIBDIR):$(DATSDIR)
VPATH=$(CLIBDIR):$(TWOCHDIR)src/:$(TWOCHDIR)include/:$(TWOCHDIR)src/OneChQS/:$(TWOCHDIR)src/TwoChQS/

##
## New pattern rule for "c" and "cpp" files
##
#%.o : %.c
#	$(CXX) -c -I./  $(CXXFLAGS) $(LDFLAGS) $< -o $@
%.o : %.c
	$(CXX) -c -g   $(CXXFLAGS) $< -o $@
%.o : %.cpp
	$(CXX) -c -g   $(CXXFLAGS) $< -o $@
##
## $^  all dependencies
## $<  just the first dependency
## $@  name of the target



##DM_NRG_SpecDens : DM_NRG_SpecDens.cpp $(OBJS) $(OBJS_DMNRG)
DM_NRG : DM_NRG.cpp $(OBJS) $(OBJS_DMNRG)
	$(CXX)  $^ $(CXXFLAGS) $(LDFLAGS) -o $@ 
##
## Don't know why it doesn't put LDFLAGS after the objects
##

Conductance  : Conductance.cpp $(OBJS) $(OBJS_CONDUCTANCE)
	$(CXX)  $^ $(CXXFLAGS) $(LDFLAGS) -o $@ 

clean :
	-rm *.o *~

