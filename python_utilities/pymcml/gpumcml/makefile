#
# Makefile for GPUMCML on Linux
#
# It generates four executables: gpumcml.sm_[11 | 12 | 13 | 20].
#
# They are optimized for GPUs with compute capability 1.1, 1.2, 1.3 and 2.0
# respectively. Note that a program can run on a GPU with compute capability
# that is higher than the one it is optimized for.
#

######################################################################
# Configurations
######################################################################

# To disable register usage constraint, set it to 0.
MAX_REG_COUNT_SM20 := 36
MAX_REG_COUNT_SM13 := 0
MAX_REG_COUNT_SM12 := 0
MAX_REG_COUNT_SM11 := 0

# To compile 32-bit binaries on a 64-bit machine, set this to 1.
COMPILE_32 := 0

######################################################################
# Compiler flags
######################################################################

CUDA_INSTALL_ROOT ?= /usr/

NVCC := $(CUDA_INSTALL_ROOT)/bin/nvcc
GCC := g++

HP_64 = $(shell uname -m | grep 64)
ifeq "$(strip $(HP_64))" ""
  # 32-bit architecture
  LIBCUDART_DIR := -L$(CUDA_INSTALL_ROOT)/lib
else
  ifeq "$(strip $(COMPILE_32))" "1"
    # 64-bit architecture, but we compile 32-bit binary.
    LIBCUDART_DIR := -L$(CUDA_INSTALL_ROOT)/lib
  else
    LIBCUDART_DIR := -L$(CUDA_INSTALL_ROOT)/lib64
  endif
endif

CFLAGS := -Wall -O3
LDFLAGS := $(LIBCUDART_DIR)

NVCC_FLAGS := -I$(CUDA_INSTALL_ROOT)/include
NVCC_FLAGS += -O3 -DUNIX -use_fast_math
NVCC_FLAGS += -Xptxas -v
#export PTX_FLAGS := $(PTXAS_FLAGS) -maxrregcount 64 -mem -v

ifeq "$(strip $(COMPILE_32))" "1"
  CFLAGS += -m32
  LDFLAGS += -m32
  NVCC_FLAGS += -m32
endif

######################################################################
# Targets
######################################################################

PROG_BASE := gpumcml
PROG_SM20 := $(PROG_BASE)

default: $(PROG_SM20)

######################################################################
# Source files
######################################################################

CU_SRCS := $(PROG_BASE)_main.cu \
		$(PROG_BASE)_kernel.cu	\
		$(PROG_BASE)_mem.cu \
		$(PROG_BASE)_rng.cu \
		$(PROG_BASE)_kernel.h	\
		$(PROG_BASE).h

######################################################################
# GPU code compilation rules for Compute Capability 2.0
######################################################################

NVCC_FLAGS_SM20 := $(NVCC_FLAGS) $(GENCODE_SM20)
ifneq "$(strip $(MAX_REG_COUNT_SM20))" "0"
NVCC_FLAGS_SM20 += -maxrregcount $(MAX_REG_COUNT_SM20)
endif
NVCC_FLAGS_SM20 += -keep # -deviceemu

$(PROG_SM20): $(PROG_SM20).cu_o gpumcml_io.o
	-$(GCC) $(LDFLAGS) -o $@ $^ -lcudart -lpthread -fPIC

$(PROG_SM20).cu_o: $(CU_SRCS)
	-$(NVCC) $(NVCC_FLAGS_SM20) -o $@ -c $<

clean_sm_20:
	-$(NVCC) $(NVCC_FLAGS_SM20) -o $(PROG_SM20).cu_o -c $(PROG_BASE)_main.cu -clean

######################################################################
# C code compilation rules
######################################################################

%.o: %.c
	$(GCC) $(CFLAGS) -c $<

%.o: %.cpp
	$(GCC) $(CFLAGS) -c $<

######################################################################

clean: clean_sm_20 clean_sm_13 clean_sm_12 clean_sm_11
	rm -f $(PROG_SM20) $(PROG_SM13) $(PROG_SM12) $(PROG_SM11) *.o *.cu_o

######################################################################

# Disable implicit rules on .cu files.
%.cu: ;

.PHONY: default clean clean_sm_20 clean_sm_13 clean_sm_12 clean_sm_11

