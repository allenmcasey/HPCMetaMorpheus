#
# Copyright (c) 2019      University of Houston. All rights reserved.
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER$
#
include ../../Makefile.defs

HEADER =  FdrAnalysisEngine.h  FdrAnalysisResults.h  FdrInfo.h


OBJ = 	  FdrAnalysisEngine.o  FdrAnalysisResults.o  FdrInfo.o


.SUFFIXES: .o.cpp

all:
	$(MAKE) $(MMORPHEUS_LIB)

FdrInfo.o: FdrInfo.cpp $(HEADER)
	$(CXX) $(CXXFLAGS) FdrInfo.cpp FdrInfo.pb.cc -c -lprotobuf -pthread

$(MMORPHEUS_LIB): $(HEADER) $(OBJ) 
	$(AR) $(AR_FLAGS) $(MMORPHEUS_LIB_DIR)/$(MMORPHEUS_LIB) $(OBJ)

$(OBJ):	$(HEADER)
.cpp.o:	$(HEADER)
	$(CXX) $(CXXFLAGS) -c $<

clean:  
	rm -rf *.o *~

proto:
	protoc -I=. --cpp_out=. FdrInfo.proto
