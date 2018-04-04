CC=mpic++
CPPFLAGS= -w -O3 -c -fpermissive
LDFLAGS= -O3
SRCS=  TaxonSet.cpp Tree.cpp Random.cpp SequenceAlignment.cpp \
	StateSpace.cpp ZippedSequenceAlignment.cpp SubMatrix.cpp \
	GTRSubMatrix.cpp linalg.cpp Chrono.cpp BranchProcess.cpp \
	GammaBranchProcess.cpp RateProcess.cpp DGamRateProcess.cpp ProfileProcess.cpp \
	ExpoConjugateGTRSBDPProfileProcess.cpp ExpoConjugateGTRPhyloProcess.cpp ExpoConjugateGTRSubstitutionProcess.cpp ExpoConjugateGTRMixtureProfileProcess.cpp \
	MixtureProfileProcess.cpp MatrixMixtureProfileProcess.cpp GTRMixtureProfileProcess.cpp MatrixSBDPProfileProcess.cpp MatrixProfileProcess.cpp SBDPProfileProcess.cpp \
	SubstitutionProcess.cpp Propagate.cpp \
	 MatrixSubstitutionProcess.cpp \
	GeneralPathSuffStatMatrixSubstitutionProcess.cpp \
	PhyloProcess.cpp MatrixPhyloProcess.cpp \
	GeneralPathSuffStatMatrixPhyloProcess.cpp GeneralPathSuffStatMatrixMixtureProfileProcess.cpp \
	RASCATGTRSBDPGammaPhyloProcess.cpp \
	Bipartition.cpp BipartitionList.cpp Consensus.cpp TaxaParameters.cpp PBTree.cpp TreeList.cpp PolyNode.cpp correl.cpp correlation.cpp \
	NNI.cpp \
	DPProfileProcess.cpp GTRProfileProcess.cpp ExpoConjugateGTRProfileProcess.cpp 



#FiniteProfileProcess.cpp PoissonSubstitutionProcess.cpp\
#MatrixDPProfileProcess.cpp MatrixFiniteProfileProcess.cpp \
#CodonSequenceAlignment.cpp CodonStateSpace.cpp 
#sure
#OneProfileProcess.cpp MatrixOneProfileProcess.cpp
#ExpoConjugateGTRFiniteProfileProcess.cpp ExpoConjugateGTRDPProfileProcess.cpp 
#PoissonMixtureProfileProcess.cpp PoissonPhyloProcess.cpp PoissonFiniteProfileProcess.cpp PoissonDPProfileProcess.cpp #PoissonSBDPProfileProcess.cpp \
#CodonSubMatrix.cpp 
#CodonMutSelProfileProcess.cpp AACodonMutSelProfileProcess.cpp 
#RASCATGammaPhyloProcess.cpp RASCATSBDPGammaPhyloProcess.cpp \
#RASCATFiniteGammaPhyloProcess.cpp RASCATGTRFiniteGammaPhyloProcess.cpp RASCATGTRDPGammaPhyloProcess.cpp \
#AACodonMutSelFinitePhyloProcess.cpp CodonMutSelFinitePhyloProcess.cpp\
#CodonMutSelSBDPPhyloProcess.cpp AACodonMutSelSBDPPhyloProcess.cpp \
	

OBJS=$(patsubst %.cpp,%.o,$(SRCS))
ALL_SRCS=$(wildcard *.cpp)
ALL_OBJS=$(patsubst %.cpp,%.o,$(ALL_SRCS))

PROGSDIR=../data
ALL= pb_mpi #readpb_mpi tracecomp bpcomp cvrep
PROGS=$(addprefix $(PROGSDIR)/, $(ALL))

.PHONY: all clean
all: $(PROGS)

# Rules for generate the dependencies automatically

%.d: %.cpp
	@echo "Generating dependencies for $<..."; \
	 set -e; rm -f $@; $(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; rm -f $@.$$$$


# Rules to create .o files from .cpp files
%.o: %.cpp %.d
	$(CC) -c $(CPPFLAGS) $*.cpp

# Include the dependencies unless the request was to clean eveything up
ifneq ($(MAKECMDGOALS),clean)
-include $(ALL_OBJS:.o=.d)
endif

# Targets

#$(PROGSDIR)/cvrep : CVrep.o $(OBJS)
#	$(CC) CVrep.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

#$(PROGSDIR)/jackknife: JackKnife.o $(OBJS)
#	$(CC) JackKnife.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

#$(PROGSDIR)/readpb_mpi: ReadPB.o $(OBJS)
#	$(CC) ReadPB.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

$(PROGSDIR)/pb_mpi: PB.o $(OBJS)
	$(CC) PB.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

#$(PROGSDIR)/tracecomp: TraceComp.o $(OBJS)
#	$(CC) TraceComp.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

#$(PROGSDIR)/bpcomp: BPCompare.o $(OBJS)
#	$(CC) BPCompare.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

clean:
	-rm -f *.o *.d *.d.*
	-rm -f $(PROGS)

