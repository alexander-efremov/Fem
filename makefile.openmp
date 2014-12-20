# hello from MAC

# Location of the CUDA Toolkit
CUDA_PATH  ?= /usr/local/cuda-6.5
NVCC := $(CUDA_PATH)/bin/nvcc -ccbin g++

# internal flags
NVCCFLAGS   := -m64 -O2
CCFLAGS     := -m64 -O2
LDFLAGS     := 


ifeq ($(flagman), 1)
   NVCCFALGS += -DFLAGMAN
endif

# Extra user flags
EXTRA_NVCCFLAGS   ?=
EXTRA_LDFLAGS     ?= -L./lib
EXTRA_CCFLAGS     ?=

# Debug build flags
ifeq ($(dbg), 1)
      NVCCFLAGS += -g -G -DDEBUG 
endif

ALL_CCFLAGS :=
ALL_CCFLAGS += $(NVCCFLAGS)
ALL_CCFLAGS += $(EXTRA_NVCCFLAGS)
ALL_CCFLAGS += $(addprefix -Xcompiler ,$(CCFLAGS))
ALL_CCFLAGS += $(addprefix -Xcompiler ,$(EXTRA_CCFLAGS))

ALL_LDFLAGS :=
ALL_LDFLAGS += $(ALL_CCFLAGS)
ALL_LDFLAGS += $(addprefix -Xlinker ,$(LDFLAGS))
ALL_LDFLAGS += $(addprefix -Xlinker ,$(EXTRA_LDFLAGS))

INCLUDES  := -I./FemBase -I./gtest-1.7.0 -I./include
LIBRARIES :=
LIBRARIES += -lgtest -lgtest_main

################################################################################

GENCODE_SM20    := -gencode arch=compute_20,code=sm_20
GENCODE_FLAGS   ?= $(GENCODE_SM20)

################################################################################

OBJ := obj

# Target rules
all: build

build: test_fixture
	
$(OBJ)/compute_density_base.o: src/LowOrdOper.cpp
	-mkdir -p $(OBJ)
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

$(OBJ)/integ_und_gor_chan.o: src/IntegUndGorChan.cpp
	-mkdir -p $(OBJ)
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

$(OBJ)/init_and_boundary_data.o: src/InitAndBoundaryData.cpp
	-mkdir -p $(OBJ)
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<
	
$(OBJ)/special_print.o: src/SpecialPrint.cpp
	-mkdir -p $(OBJ)
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

$(OBJ)/space_volume.o: src/SpaceVolume.cpp
	-mkdir -p $(OBJ)
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

$(OBJ)/integ_und_rig_ang_tr.o: src/integUndRigAngTr.cpp
	-mkdir -p $(OBJ)
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

$(OBJ)/one_cell_integ.o: src/OneCellInteg.cpp
	-mkdir -p $(OBJ)
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<
	
$(OBJ)/origin_print.o: src/OriginPrint.cpp
	-mkdir -p $(OBJ)
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

$(OBJ)/compute_density.o: src/compute_density.cpp
	-mkdir -p $(OBJ)
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

$(OBJ)/tests.o: src/tests.cpp
	-mkdir -p $(OBJ)
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

$(OBJ)/gtest_main.o: src/gtest_main.cc
	-mkdir -p $(OBJ)
	$(NVCC) $(INCLUDES) $(ALL_CCFLAGS) $(GENCODE_FLAGS) -o $@ -c $<

test_fixture: $(OBJ)/special_print.o $(OBJ)/origin_print.o $(OBJ)/one_cell_integ.o $(OBJ)/space_volume.o $(OBJ)/integ_und_rig_ang_tr.o $(OBJ)/gtest_main.o $(OBJ)/integ_und_gor_chan.o $(OBJ)/compute_density_base.o $(OBJ)/init_and_boundary_data.o $(OBJ)/compute_density.o $(OBJ)/tests.o
	$(NVCC) $(ALL_LDFLAGS) $(GENCODE_FLAGS) -o $@ $+ $(LIBRARIES)

run: build
	./test_fixture

clean:
	rm -f test_fixture $(OBJ)/*.o

clobber: clean
