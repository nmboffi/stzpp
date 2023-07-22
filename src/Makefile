# Manually set the compiler
CXX = mpicxx -Wno-long-long -fopenmp

# Compiler flags
CXXFLAGS = -Wall -lm -O3 -pedantic -std=c++11 -L../mg_3d_template/ -fPIC

# Object files and corresponding source files
OBJS = common.o shear_sim_3d.o trans_sim_3d.o \
	   ../mg_3d_template/buffer.o ../mg_3d_template/common.o ../mg_3d_template/geometry.o

SHEAR_SIM_OBJS     = shear_sim.o
TRANS_SIM_OBJS     = trans_shear_sim.o
SC_OBJS            = shear_compare.o
SCT_OBJS           = shear_compare_trans.o

SRCS = $(OBJS:%.o=%.cpp) shear_sim.cpp trans_shear_sim.cpp shear_compare.cpp shear_compare_trans.cpp mat3.cc

# List of all executable files to compile
EXECS = shear_sim trans_shear_sim shear_compare shear_compare_trans

# How to create all of the executables
all: 
	@echo Making executables: $(EXECS)
	$(MAKE) executables
   
# Target listing all the executables to be built
executables: $(OBJS) $(SHEAR_SIM_OBJS) $(TRANS_SIM_OBJS) $(EXECS)

# Infer dependencies and create specific dependency files
# See bottom of http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/o
# Put dependencies in .d
DEPDIR = .d

# Create DEPDIR if it does not exist, creating all needed directories recursively,
# and toss the output in the bitbucket
$(shell mkdir -p $(DEPDIR) >/dev/null)

# Explanation as follows:
# -MT $@: put the name of the target in the generated dependency file
# -MMD: generate dependency information as a side-effect of compilation, not instead,
#  and also omit system headers
# -MP: Add target for each prerequisite in the list
# -MF $(DEPDIR)/$*.Td: first create the dep files as temporary files with a .Td extension
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td

# Compile using CXX, with the DEPFLAGS and CXXFLAGS, for the TARGET_ARCH architecture
COMPILE.cpp = $(CXX) $(DEPFLAGS) $(CXXFLAGS) $(TARGET_ARCH) -c

# After finishing compilation, move the temporary files over to "real" dependency files
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d

# Overwrite the original rule for creating object files
%.o : %.cpp

# For every object file, the dependency file is a prerequisite.
# This guarantees that if a .d file has been updated, the .o file will be recompiled,
# and if the .hh file has changed, the .d file will change too.
# Then just compile all the files with any OUTPUT_OPTIONS (none) and output it as the name of the target
%.o : %.cpp $(DEPDIR)/%.d
	$(COMPILE.cpp) $<
	$(POSTCOMPILE)

# Empty rule for each dependency file so make won't fail if one doesn't exist
$(DEPDIR)/%.d: ;

# Declare each dependency file as PRECIOUS so they wont be deleted as intermediate files
.PRECIOUS: $(DEPDIR)/%.d

# Include the dependency file for every source file
-include $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS)))

# Compile the shear_sim executable
shear_sim: $(OBJS) $(SHEAR_SIM_OBJS)
	@echo Making shear_sim
	$(CXX) $^ -o $@ $(CXXFLAGS) -llapack -lblas

# Compile the trans_sim executable
trans_shear_sim: $(OBJS) $(TRANS_SIM_OBJS)
	@echo Making trans_sim 
	$(CXX) $^ -o $@ $(CXXFLAGS) -llapack -lblas

# Compile the shear_compare executable
shear_compare: $(OBJS) $(SC_OBJS)
	@echo Making shear_compare 
	$(CXX) $^ -o $@ $(CXXFLAGS) -llapack -lblas

# Compile the shear_compare executable
shear_compare_trans: $(OBJS) $(SCT_OBJS)
	@echo Making shear_compare_trans 
	$(CXX) $^ -o $@ $(CXXFLAGS) -llapack -lblas


.PHONY: clean cleanobj cleanexecs cleandat all executables depend

# How to remove old files
cleanobj:
	rm *.o

cleandep:
	rm -r .d

cleandat:
	rm *.dat

cleanexecs:
	rm $(EXECS)

clean: cleanobj cleanexecs cleandep
