# COMPILER AND LINKER OPTIONS
CC = g++-5
CPP = g++-5

NVCC = /usr/local/cuda/bin/nvcc
COMPUTE_CAPABILITY = 61
NVCCFLAGS  = -gencode arch=compute_$(COMPUTE_CAPABILITY),code=sm_$(COMPUTE_CAPABILITY)
NVCCFLAGS += -m64 -std=c++11  -ccbin $(CPP) 
COMPILER_OPTIONS = -Wall -Wno-deprecated-declarations -pthread
NVCCFLAGS += -g -G -Xcompiler "$(COMPILER_OPTIONS)" -O2

INCLUDES = -I/usr/local/cuda/samples/common/inc

LDFLAGS = -lglut -lGLU -lGL -lGLEW -lm -lpthread
LDFLAGS += -L/usr/local/cuda/lib64

# PATHS CONFIGURATION
CIBLE = main
SRC_DIR = ./src
OBJ_DIR = ./build

CPP_FILES = $(wildcard $(SRC_DIR)/*.cpp)
CU_FILES  = $(wildcard $(SRC_DIR)/*.cu)

OBJS  = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(CPP_FILES))
OBJS += $(patsubst $(SRC_DIR)/%.cu,$(OBJ_DIR)/%.cu.o,$(CU_FILES))

# AUTOMATICALLY HANDLE HEADER INCLUSION
$(OBJ_DIR)/%.d : $(SRC_DIR)/%.cpp 
	$(CC) $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -MM -MF $@ -MP $<
$(OBJ_DIR)/%.cu.d : $(SRC_DIR)/%.cu
	$(CC) $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -MM -MF $@ -MP $<

#BUILDING RULES
.PHONY: default
default: $(CIBLE)

$(CIBLE): $(OBJS)
	$(NVCC) $^ -o $@ $(LDFLAGS) -arch=sm_$(COMPUTE_CAPABILITY)
$(OBJ_DIR)/%.cu.o : $(SRC_DIR)/%.cu $(OBJ_DIR)/%.cu.d
	@mkdir -p $(OBJ_DIR)
	$(NVCC) $(NVCCFLAGS) $(INCLUDES) -dc -o $@ $<
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp $(OBJ_DIR)/%.d
	@mkdir -p $(OBJ_DIR)
	$(NVCC) $(NVCCFLAGS) $(INCLUDES) -dc -o $@ $<

clean:
	rm -f  $(CIBLE) build/*