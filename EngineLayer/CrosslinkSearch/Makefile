#
# Copyright (c) 2019      University of Houston. All rights reserved.
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER$
#
include ../../Makefile.defs

HEADER =  BestPeptideScoreNotch.h  \
          Crosslinker.h            \
          CrosslinkSpectralMatch.h \
          CrosslinkedPeptides.h    \
          CrosslinkSearchEngine.h  \
          PsmCrossType.h


OBJ =     BestPeptideScoreNotch.o  \
          Crosslinker.o            \
          CrosslinkSpectralMatch.o \
          CrosslinkedPeptides.o    \
          CrosslinkSearchEngine.o  \
          PsmCrossType.o



.SUFFIXES: .o.cpp

all:
	$(MAKE) $(MMORPHEUS_LIB)

$(MMORPHEUS_LIB): $(HEADER) $(OBJ)
	$(AR) $(AR_FLAGS) $(MMORPHEUS_LIB_DIR)/$(MMORPHEUS_LIB) $(OBJ)

$(OBJ): $(HEADER)
.cpp.o: $(HEADER)
	$(CXX) $(CXXFLAGS) -I../FdrAnalysis -I../../mzlib/Proteomics/ProteolyticDigestion -I../../mzlib/Proteomics/Fragmentation -c $< 

clean:
	rm -rf *.o *~

flatbuffers:
	flatc --cpp -I ../FdrAnalysis -I ../../mzlib/Proteomics/ProteolyticDigestion -I ../../mzlib/Proteomics/Fragmentation CrosslinkSpectralMatch.fbs
