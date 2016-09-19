# Makefile for OS-X

MAYA_LOCATION = /Applications/Autodesk/maya

PROJ1 = cageDeformer
PROJ2 = cageDeformerARAP
SOURCE1 = ./$(PROJ1)/$(PROJ1).cpp
SOURCE2 = ./$(PROJ2)/$(PROJ2).cpp
OBJECT1 = ./$(PROJ1).o
OBJECT2 = ./$(PROJ2).o
ASM1 = ./$(PROJ1).s
ASM2 = ./$(PROJ2).s

INCLUDES = -I$(MAYA_LOCATION)/devkit/include/ -I./ -I../ -I/usr/local/include/eigen3/
LIBS = -L$(MAYA_LOCATION)/Maya.app/Contents/MacOS -lOpenMaya -lOpenMayaAnim -lOpenMayaRender -lOpenMayaUI -lFoundation
ISYSROOT = /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk

#compiler options
CC = /usr/bin/g++
LD = /usr/bin/g++
CLANG = /usr/bin/clang

CFLAGS = -DCC_GNU_ -DOSMac_ -DOSMacOSX_  \
	-DOSMac_MachO_ -D_LANGUAGE_C_PLUS_PLUS \
	-include "$(MAYA_LOCATION)/devkit/include/maya/OpenMayaMac.h" \
	-DMAC_PLUGIN -D_BOOL -DREQUIRE_IOSTREAM -O2 -mavx
#CFLAGS +=  -fopenmp
C++FLAGS = $(CFLAGS) $(WARNFLAGS) $(ERROR_FLAGS) -fno-gnu-keywords
LDFLAGS = $(C++FLAGS)
DYNLIB_LOCATION	= $(MAYA_LOCATION)/Maya.app/Contents/MacOS

LDFLAGS	+= -isysroot $(ISYSROOT) $(ARCH_FLAGS) -headerpad_max_install_names -bundle

LDFLAGS += -Wl,-exported_symbol,__Z16initializePlugin7MObject \
          -Wl,-exported_symbol,__Z18uninitializePlugin7MObject

LREMAP = -Wl,-executable_path,"$(DYNLIB_LOCATION)"
LDFLAGS += -L"$(DYNLIB_LOCATION)" $(LREMAP)

.PHONY: all install clean


all: $(PROJ1) $(PROJ2)

$(PROJ1): $(OBJECT1)
	        $(LD) $(LDFLAGS) $(LIBS) $^ -o $@.bundle

$(PROJ2): $(OBJECT2)
	        $(LD) $(LDFLAGS) $(LIBS) $^ -o $@.bundle

$(OBJECT1): $(ASM1)
	        $(CLANG) -c $^

$(OBJECT2): $(ASM2)
	        $(CLANG) -c $^

$(ASM1): $(SOURCE1)
	        $(CC) -S $(INCLUDES) $(C++FLAGS) $^

$(ASM2): $(SOURCE2)
	        $(CC) -S $(INCLUDES) $(C++FLAGS) $^

install:
	        mv $(PROJ1).bundle $(PROJ2).bundle /Users/Shared/Autodesk/maya/plug-ins/

clean:
		rm -f $(OBJECT1) $(OBJECT2) $(ASM1) $(ASM2)
