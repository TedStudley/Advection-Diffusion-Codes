# C/C++ compiler
CC					:=	g++
# C/C++ compiler flags
CFLAGS			+=	-Wall -c -Iinclude -fopenmp -std=c++11

# C/C++ linker
LD					:= 	g++
# C/C++ linker flags
LDFLAGS			+=


# Source modules
MODULES			:=  
# Build directory and target
BUILD_DIR		:= 	build
EXECUTABLE	:= 	advection-1d
TARGET			:= 	$(addprefix $(BUILD_DIR)/,$(EXECUTABLE))
# Source directory and files
SRC_DIR			:=	source
SRC_DIRS		:=	$(addprefix $(SRC_DIR)/,$(MODULES)) $(SRC_DIR)
SOURCES			:=	$(foreach sdir,$(SRC_DIRS),$(wildcard $(sdir)/*.cpp))
# Object directory and files
OBJ_DIR			:= $(addprefix $(BUILD_DIR)/,obj)
OBJ_DIRS		:= $(addprefix $(OBJ_DIR)/,$(MODULES))	$(OBJ_DIR)
OBJECTS			:= $(foreach src,$(SOURCES),$(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(src)))
# Header directory and files
HDR_DIR			:= include
HDR_DIRS		:= $(addprefix $(HDR_DIR)/,$(MODULES)) $(HDR_DIR)
HEADERS			:= $(foreach hdir,$(HDR_DIRS),$(wildcard $(hdir)/*.h))

vpath %.cpp	$(patsubst ' ',':',$(SRC_DIRS))
vpath %.o	$(patsubst ' ',':',$(OBJ_DIRS))

define make-goal
$1/%.o %.cpp
	$(CC) $(CFLAGS) -c $$< -o $$@
endef

.PHONY: all checkdirs clean style

all: checkdirs $(TARGET)

$(OBJ_DIR)/%.o : %.cpp
	@echo =====\($(EXECUTABLE)\)====================== Building $<
	@$(CC) $(CFLAGS) -c $< -o $@
	
$(TARGET) : $(OBJECTS)
	@echo =====\($(EXECUTABLE)\)====================== Linking $@
	@$(LD) $(LDFLAGS) $^ -o $@
	
checkdirs: $(BUILD_DIR) $(OBJ_DIRS)

$(BUILD_DIR) $(OBJ_DIRS):
	@mkdir -p $@

clean:
	@rm -rf $(BUILD_DIR)
