include make.defs.oneapi

BIN = bin
OBJ_DIR = obj
RED = \033[1;31m
GREEN = \033[1;32m
BLUE = \033[1;34m
YELLOW = \033[1;33m
NC = \033[1;0m
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
LIB = -framework Accelerate
TARGET = $(BIN)/asfem_macos.exe
endif
ifeq ($(UNAME), Linux)
LIB = -llapack 
TARGET = $(BIN)/asfem_linux.exe
endif
CC = c++

#************ CLASS DIR ***********
STATIC = ../Static
KINECT = ../FP_FEM_Matrices/Kinect
OVERLAP = ../FP_FEM_Matrices/Overlap
EXTPOT = ../FP_FEM_Matrices/ExternalPot
NUMERICFEM = ../FP_FEM_Matrices/Numeric
FEM = ../FP_Model
FEM_MAT = ../FP_FEM_Matrices
ATOMIC = ../Atomic_Computations
ES = ../ElectronicStructure
TOOLS = ../Tools
GRID = ../Grids
MATRIX = ../Matrix
SYCL = ../SYCL
# ********** SOURCE *************
SOURCE += main.cpp
SOURCE += $(wildcard $(FEM_MAT)/*.cpp)
SOURCE += $(wildcard $(OVERLAP)/*.cpp)
SOURCE += $(wildcard $(KINECT)/*.cpp)
SOURCE += $(wildcard $(EXTPOT)/*.cpp)
SOURCE += $(wildcard $(NUMERICFEM)/*.cpp)
SOURCE += $(STATIC)/grid_tools.cpp
SOURCE += $(FEM)/FemFixedPoints.cpp
SOURCE += $(ATOMIC)/atomicStructure.cpp
SOURCE += $(ES)/closedShell.cpp
SOURCE += $(wildcard $(TOOLS)/*.cpp)

#****** HEADERS ************
INC += $(wildcard $(FEM_MAT)/*.hpp)
INC += $(wildcard $(OVERLAP)/*.hpp)
INC += $(wildcard $(KINECT)/*.hpp)
INC += $(wildcard $(EXTPOT)/*.hpp)
INC += $(wildcard $(NUMERICFEM)/*.hpp)
INC += $(STATIC)/grid_tools.hpp
INC += $(FEM)/FemFixedPoints.hpp
INC += $(ATOMIC)/atomicStructure.hpp
INC += $(ES)/closedShell.hpp
INC += $(wildcard $(TOOLS)/*.hpp)
INC += $(wildcard $(GRID)/*.hpp)
INC += $(wildcard $(MATRIX)/*.hpp)
#**********
VPATH = $(FEM) $(STATIC) $(FEM_MAT) $(ATOMIC) $(ES) $(TOOLS) $(OVERLAP) $(KINECT) $(EXTPOT) $(NUMERICFEM)
OBJECT = $(patsubst %,$(OBJ_DIR)/%, $(notdir $(SOURCE:.cpp=.o)))
#OBJ = $(SOURCE:.cpp=.o)
$(TARGET) : $(OBJECT)
	@echo -e "$(RED)LINKING...$(NC)"
	$(CC) -o  $@ $^ -lm $(LIB)
	@echo -e "$(YELLOW)FINISHED!!!$(NC)"
	@echo -e "$(YELLOW)SUCCESFULL!!!$(NC)"
	@echo -e "$(BLUE)OPERATING SISTEM: $(UNAME)!!!...$(NC)"
$(OBJ_DIR)/%.o : %.cpp
	@echo -e "$(GREEN)COMPILING...$(NC)"
	@echo -e  "\033[0;32m [OK] \033[0m       \033[0;33m Compiling:\033[0m" $<
	$(CC) $(CFLAGS) -c  $< -o $@
print-%  : ; @echo $* = $($*)


.PHONY : help clean

clean :
	rm -f $(OBJECT) $(BIN)/$(TARGET)
