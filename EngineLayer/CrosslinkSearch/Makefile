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


OBJ = 	  BestPeptideScoreNotch.o  \
	  Crosslinker.o            \
	  CrosslinkSpectralMatch.pb.o \
	  CrosslinkSpectralMatch.o \
	  CrosslinkedPeptides.o    \
	  CrosslinkSearchEngine.o  \
	  PsmCrossType.o



.SUFFIXES: .o.cpp

all:
	$(MAKE) $(MMORPHEUS_LIB)

#CrosslinkSpectralMatch.o: CrosslinkSpectralMatch.cpp $(HEADER)
#	$(CXX) $(CXXFLAGS) -I/../FdrAnalysis -I../../mzlib/Proteomics/ProteolyticDigestion -I../../mzlib/Proteomics/Fragmentation CrosslinkSpectralMatch.cpp CrosslinkSpectralMatch.pb.cc -c -lprotobuf

$(MMORPHEUS_LIB): $(HEADER) $(OBJ)
	$(AR) $(AR_FLAGS) $(MMORPHEUS_LIB_DIR)/$(MMORPHEUS_LIB) $(OBJ) -lprotobuf -pthread

$(OBJ):	$(HEADER)
.cpp.o:	$(HEADER)
	$(CXX) $(CXXFLAGS) -I../FdrAnalysis -I../../mzlib/Proteomics/ProteolyticDigestion -I../../mzlib/Proteomics/Fragmentation -c $< -lprotobuf -pthread

.cc.o: $(HEADER)
	$(CXX) $(CXXFLAGS) -I../FdrAnalysis -I../../mzlib/Proteomics/ProteolyticDigestion -I../../mzlib/Proteomics/Fragmentation -c $< -lprotobuf -pthread

clean:  
	rm -rf *.o *~

proto:
	protoc -I=. --cpp_out=. --proto_path=../FdrAnalysis --proto_path=../../mzlib/Proteomics/ProteolyticDigestion --proto_path=../../mzlib/Proteomics/Fragmentation CrosslinkSpectralMatch.proto
