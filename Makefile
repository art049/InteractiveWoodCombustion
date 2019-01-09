# COMPILER AND LINKER OPTIONS
CC = g++-5
CPP = g++-5
CPPFLAGS = -std=c++11
NVCC = /usr/local/cuda/bin/nvcc
COMPUTE_CAPABILITY = 61
#COMPUTE_CAPABILITY = 50

NVCCFLAGS  = -gencode arch=compute_$(COMPUTE_CAPABILITY),code=sm_$(COMPUTE_CAPABILITY)
NVCCFLAGS += -m64 -std=c++11  -ccbin $(CPP) 
COMPILER_OPTIONS = -Wall -Wno-deprecated-declarations -pthread -DGL_GLEXT_PROTOTYPES
NVCCFLAGS += -g -G -Xcompiler "$(COMPILER_OPTIONS)" -O2

INCLUDES = -I/usr/local/cuda/samples/common/inc 

LDFLAGS = -lglut -lGLU -lGL -lGLEW -lm -lpthread
LDFLAGS += -L/usr/local/cuda/lib64

# PATHS CONFIGURATION
CIBLE = main
SRC_DIR = ./src
OBJ_DIR = ./build

CPP_FILES = $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/**/*.cpp)
CU_FILES  = $(wildcard $(SRC_DIR)/*.cu) $(wildcard $(SRC_DIR)/**/*.cu)

OBJS  = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(CPP_FILES))
OBJS += $(patsubst $(SRC_DIR)/%.cu,$(OBJ_DIR)/%.cu.o,$(CU_FILES))

# AUTOMATICALLY HANDLE HEADER INCLUSION
$(OBJ_DIR)/%.d : $(SRC_DIR)/%.cpp 
	@mkdir -p $(OBJ_DIR)/physics $(OBJ_DIR)/cuda_common
	$(CC) $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -MM -MF $@ -MP $<
$(OBJ_DIR)/%.cu.d : $(SRC_DIR)/%.cu
	@mkdir -p $(OBJ_DIR)/physics $(OBJ_DIR)/cuda_common
	$(NVCC) $(NVCCFLAGS) $(TARGET_ARCH) -M $< > $@
-include $($(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.d,$(CPP_FILES)))
-include $($(patsubst $(SRC_DIR)/%.cu,$(OBJ_DIR)/%.cu.d,$(CU_FILES)))

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
	rm -rf  $(CIBLE) build/* 