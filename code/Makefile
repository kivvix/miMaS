# miMaS compilation
PROJET = miMaS

DEBUG_MODE ?= N

include Makefile.config

SRC_dir    = source
HEADER_dir = include/$(PROJET)
OBJ_dir    = obj

INC_BOOST = /usr/local/Cellar/boost/1.71.0/include
INC_FFTW  = /usr/local/Cellar/fftw/3.3.8_1/include
LIB_FFTW  = fftw3

INCLUDES = $(INC_BOOST) $(INC_FFTW) $(HEADER_dir)/..
LIBS = $(PROJET) $(LIB_FFTW) m

### Lists of files and flags for compilation ##########################
SRC = $(wildcard $(SRC_dir)/*.cc)
OBJ = $(SRC:$(SRC_dir)/%.cc=$(OBJ_dir)/%.o)
INC = $(wildcard $(HEADER_dir)/*.h)

# compilation flags
GLLIBS = -lm

ifeq ($(DEBUG_MODE),Y)
	CFLAGS += -O0 -DDEBUG_MODE
	DEBUGFLAGS = -p
	
	COMPIL_mode = "\033[1mCompilation en mode debug\033[0m"
else
	CFLAGS += -O3

	COMPIL_mode = "\033[1mCompilation en mode release\033[0m"
endif

#linker
LD = $(CXX) -fPIC -shared

LDFLAGS = $(CFLAGS)

### BUILD RULES #######################################################

all : $(LIB_dir)/lib$(PROJET).so
	@echo "all your base are belong to us"
	@echo "Library build : \033[90mlib$(PROJET).so\033[0m"

$(LIB_dir)/lib$(PROJET).so : $(LIB_dir) $(OBJ_dir) $(OBJ)
	@echo $(COMPIL_mode)
	@echo "\033[36m$(notdir $@)\033[0m"
	$(LD) $(OBJ) $(LDFLAGS) $(GLLIBS) -o $@
	@echo "\033[90mCompilation finie.\033[0m"

$(OBJ_dir)/%.o : $(SRC_dir)/%.cc $(HEADER_dir)/%.h
	@echo "\033[95m$@\033[0m"
	$(CXX) $(CFLAGS) -o $@ -c $<

$(OBJ_dir) :
	@echo "\033[92mmake $@ directory\033[0m"
	mkdir -p $@
$(LIB_dir) :
	@echo "\033[92mmake $@ directory\033[0m"
	mkdir -p $@

### cleanning rules

clean :
	@echo "\033[41;97m ** Suppression des fichiers objets ** \033[0m"
	rm -f $(OBJ_dir)/*.o

mrproper : clean
	@echo "\033[31;1m Suppression de la bibliothèque \033[0m"
	rm -f $(LIB_dir)/lib$(PROJET).so

nuke : mrproper
	@$(echo) "  ,-*\n (_)\n"
	rm -f *.out

coffee :
	@echo " (\n  )\nc[]"

### compilation rules

%.out : %.cc $(LIB_dir)/lib$(PROJET).so $(INC)
	@echo "\033[95m$@\033[0m"
	$(CXX) $(CFLAGS) $(DEBUGFLAGS) -g $< $(foreach inc,$(INCLUDES),-I $(inc)) -L $(LIB_dir) $(foreach lib,$(LIBS),-l$(lib)) -o $@
	@echo "\033[1mcompilation done\033[0m"


