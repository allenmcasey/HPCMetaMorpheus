#
# Copyright (c) 2019      University of Houston. All rights reserved.
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER$
#
include ../../Makefile.defs

HEADER =  XLSearchTask.h  \
	  XLSearchParameters.h

OBJ = 	 XLSearchTask.o  \
	 XLSearchParameters.o



.SUFFIXES: .o.cpp

all:    
	$(MAKE) $(MMORPHEUS_LIB) 


$(MMORPHEUS_LIB): $(HEADER) $(OBJ) 
	$(AR) $(AR_FLAGS) $(MMORPHEUS_LIB_DIR)/$(MMORPHEUS_LIB) $(OBJ)

$(OBJ):	$(HEADER)
.cpp.o:	$(HEADER)
	$(CXX) $(CXXFLAGS) -I../../EngineLayer/FdrAnalysis -I../../mzlib/Proteomics/ProteolyticDigestion -I../../mzlib/Proteomics/Fragmentation -c $<

.cc.o: $(HEADER)
	$(CXX) $(CXXFLAGS) -I../../EngineLayer/FdrAnalysis -I../../mzlib/Proteomics/ProteolyticDigestion -I../../mzlib/Proteomics/Fragmentation -c $<

clean:  
	rm -rf *.o *~
