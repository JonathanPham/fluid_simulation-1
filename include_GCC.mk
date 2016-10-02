CC  = gcc
CXX = g++
FC  = gfortran
LINKER = $(CXX)

ANSI_CFLAGS  = -ansi
ANSI_CFLAGS += -std=c++0x
ANSI_CFLAGS += -pedantic
ANSI_CFLAGS += -Wextra
ANSI_CFLAGS += -fstrict-aliasing
ANSI_CFLAGS += -Wunsafe-loop-optimizations

CFLAGS   = -O3 -Wno-format -fopenmp -Wall $(ANSI_CFLAGS)
CXXFLAGS = $(CFLAGS)
FCFLAGS  = 
CPPFLAGS = -std=c++0x
LFLAGS   =
DEFINES  = -D_GNU_SOURCE -DNDEBUG
INCLUDES =
LIBS     = -fopenmp
