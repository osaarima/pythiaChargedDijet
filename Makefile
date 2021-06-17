PROGRAM       = pythiaChargedDijet

version       = JTKT
CXX           = g++
#CXXFLAGS      = -O -Wall -g -Wno-deprecated -bind_at_load -D$(version)
CXXFLAGS      = -O -Wall -g -Wno-deprecated -Wno-misleading-indentation -D$(version) #-ggdb
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
#############################################
# -bind_at_load helps to remove linker error
############################################
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS  = $(shell root-config --libs)
CXXFLAGS += $(shell $(FASTJET)/bin/fastjet-config --cxxflags )
#LDFLAGS += -Wl,-rpath,/cvmfs/alice.cern.ch/x86_64-2.6-gnu-4.1.2/Packages/fastjet/v3.2.1_1.024-alice1-3/lib -lm  -L/cvmfs/alice.cern.ch/x86_64-2.6-gnu-4.1.2/Packages/fastjet/v3.2.1_1.024-alice1-3/lib -lfastjettools -lfastjet 
#LDFLAGS += -Wl,-rpath,/cvmfs/alice.cern.ch/x86_64-2.6-gnu-4.1.2/Packages/fastjet/v3.0.6_1.012-7/lib -lm  -L/cvmfs/alice.cern.ch/x86_64-2.6-gnu-4.1.2/Packages/fastjet/v3.0.6_1.012-7/lib -lfastjettools -lfastjet
LDFLAGS += $(shell $(FASTJET)/bin/fastjet-config --libs --plugins ) 

# When running pythia with lhapdf:
#LDFLAGS += -Wl,-rpath,$(PYTHIA8)/lib -lpythia8lhapdf5

#Note: not needed in pythia83**
#LDFLAGS += -L$(PYTHIA8)/lib -lpythia8 
INCS    += -I$(PYTHIA8)/include
#Instead:
LDFLAGS += $(shell $(PYTHIA8)/bin/pythia8-config --cppflags --libs)

#LDFLAGS += -L$(HEPMC)/lib -lHepMC 
#LDFLAGS += -Wl,-rpath -Wl,$(HEPMC)/lib
#LDFLAGS += -L$(HEPPDT)/lib -lHepPDT -lHepPID
#LDFLAGS += -Wl,-rpath -Wl,$(HEPPDT)/lib
#INCS    += -I$(HEPMC)/include
#INCS    += -I$(HEPPDT)/include
LDFLAGS += -L/home/alidock/.sw/slc7_x86-64/cgal/latest-aliroot6-user-next-root6/lib/ -lCGAL 
LDFLAGS += -L/home/alidock/.sw/slc7_x86-64/GMP/latest-aliroot6-user-next-root6/lib/ -lgmp
CXXFLAGS  += $(INCS)
LDFLAGS += $L -ldl

#HDRSDICT = src/AliJBaseCard.h src/AliJCard.h src/JHistos.h src/AliJBaseTrack.h #For pythiaChargedJet and others
HDRSDICT = src/AliJCDijetHistos.h \
           src/AliJCDijetAna.h \
           src/AliJHistogramInterface.h \
           src/AliJHistManager.h \
           src/AliJBaseTrack.h \
           src/AliJBaseCard.h \
           src/AliJCard.h \
           src/AliJPhoton.h #For pythiaChargedDijet
           
HDRS	+= $(HDRSDICT)  nanoDict.h


SRCS = $(HDRS:.h=.cxx)
OBJS = $(HDRS:.h=.o)

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS) src/AliJConst.h $(PROGRAM).C
		@echo "Linking $(PROGRAM) ..."
		$(CXX)  -lPhysics -L$(PWD) $(PROGRAM).C $(CXXFLAGS) $(OBJS) $(LDFLAGS) -o $(PROGRAM) 
		@echo "finally done"

%.cxx:

%: %.cxx
#  commands to execute (built-in):
	$(LINK.cc) $^ $(CXXFLAGS) $(LOADLIBES) $(LDLIBS) -o $@

%.o: %.cxx %.h
#  commands to execute (built-in):
	$(COMPILE.cc) $(OUTPUT_OPTION) $<


clean:
		rm -f $(OBJS) core *Dict* $(PROGRAM).o *.d $(PROGRAM) $(PROGRAM).sl

cl:  clean $(PROGRAM)

nanoDict.cc: $(HDRSDICT)
		@echo "Generating dictionary ..."
		@rm -f nanoDict.cc nanoDict.hh nanoDict.h
		@rootcint nanoDict.cc -c -D$(version) $(HDRSDICT)
