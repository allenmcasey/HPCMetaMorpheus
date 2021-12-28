#
# Copyright (c) 2019-2021    University of Houston. All rights reserved.
# $COPYRIGHT$
#
# Additional copyrights may follow
#
# $HEADER$
#
include Makefile.defs


all:    
	cd mzlib ; $(MAKE) ; cd ..
	cd EngineLayer ; $(MAKE) ; cd ..
	cd TaskLayer ; $(MAKE) ; cd ..
	cd CMD ; $(MAKE) ; cd ..

flatbuffers:
	cd mzlib/Proteomics/ProteolyticDigestion ; $(MAKE) flatbuffers ; cd ../../..
	cd mzlib/Proteomics/Fragmentation ; $(MAKE) flatbuffers ; cd ../../..
	cd EngineLayer/FdrAnalysis ; $(MAKE) flatbuffers ; cd ../..
	cd EngineLayer/CrosslinkSearch ; $(MAKE) flatbuffers ; cd ../..

clean:  
	cd mzlib ; $(MAKE) clean ; cd ..
	cd EngineLayer ; $(MAKE) clean ; cd ..
	cd TaskLayer ; $(MAKE) clean ; cd ..
	cd CMD ; $(MAKE) clean ; cd ..
	cd Test ; $(MAKE) clean ; cd ..
	cd $(MMORPHEUS_LIB_DIR) ; rm -f $(MMORPHEUS_LIB) ; cd ..
	cd $(MMORPHEUS_BIN_DIR) ; rm -f HPCMetaMorpheus ; cd ..	
	rm -rf *.o *~

links:
	@test -s $(MMORPHEUS_BIN_DIR)/elements.dat || ln -s $(MMORPHEUS_DIR)/Test/Data/elements.dat $(MMORPHEUS_BIN_DIR)
	@test -s $(MMORPHEUS_BIN_DIR)/proteases.tsv || ln -s $(MMORPHEUS_DIR)/Test/proteases.tsv $(MMORPHEUS_BIN_DIR)
	@test -s $(MMORPHEUS_BIN_DIR)/settings.toml || ln -s $(MMORPHEUS_DIR)/Test/settings.toml $(MMORPHEUS_BIN_DIR)
	@test -s "$(MMORPHEUS_BIN_DIR)/Data" || ln -s $(MMORPHEUS_DIR)/Test/Data/ $(MMORPHEUS_BIN_DIR)
	@test -s "$(MMORPHEUS_BIN_DIR)/Mods" || ln -s $(MMORPHEUS_DIR)/Test/Mods/ $(MMORPHEUS_BIN_DIR)
	@test -s "$(MMORPHEUS_BIN_DIR)/Contaminants" || ln -s $(MMORPHEUS_DIR)/Test/Contaminants/ $(MMORPHEUS_BIN_DIR)
