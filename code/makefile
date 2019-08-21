# To use this makefile, set FFTWHOME as an environment variable containing the
# path to the FFTW3 software. The contents of this directory should be the
# usual "bin include lib share". In the Comet cluster, FFTWHOME is
# automatically set if you have "module load fftw/3.3.4" in ~/.bashrc

#CC = g++
CC = mpic++
DFLAGS = -DPAR
CFLAGS =  -O3 -march=native -w -Wall
INCLUDE = -I${FFTWHOME}/include
LDFLAGS = -lm -lfftw3_mpi -lfftw3 -L${FFTWHOME}/lib
FLAGS     = -I${FFTWHOME}/include  -O3 -Wno-unused-result -Wno-write-strings
LIBS      = -lm -lfftw3_mpi -lfftw3 -O3  -L${FFTWHOME}/lib
DLIBS      = -lm -lfftw3_mpi -lfftw3  -L${FFTWHOME}/lib
DFLAGS     = -I${FFTWHOME}/include -Wno-unused-result -Wno-write-strings


SRCS = main.cpp calc_poly_density.cpp homopolymer_discrete.cpp \
       io_utils.cpp integrate_utils.cpp ran2.cpp array_utils.cpp calc_h.cpp \
       fft_wrappers.cpp calc_debye.cpp initialize.cpp \
       read_input.cpp diblock_discrete.cpp 1s_update.cpp \
       Euler_update.cpp simulate.cpp sphere_utils.cpp \
       grafted_nanoparticle.cpp channel.cpp

DEBUGFOLDER := objects/debug
MAINOBJFOLDER := objects/main
DIM       = $(shell grep -e "^\#define Dim" globals.h | awk '{print $$3}')
# EXE := $(shell git describe --always)-$(DIM)d.out
EXE := a_ft-gnp
TAGS := .tags


OBJS = ${SRCS:.cpp=.o}
DOBJS = $(addprefix $(DEBUGFOLDER)/, ${OBJS} )
OBJS_MAKE = $(addprefix $(MAINOBJFOLDER)/, ${OBJS} )

$(DEBUGFOLDER)/%.o: %.cpp
	        @mkdir -p $(dir $@) 
		${CC} -g ${DFLAGS} -c $< -o $@ -DPAR

$(MAINOBJFOLDER)/%.o: %.cpp
	        @mkdir -p $(dir $@) 
	        ${CC} ${CFLAGS} ${DFLAGS} -c $< -o $@  -DPAR

all: $(EXE) 

.cpp.o:
	${CC} ${CFLAGS} ${DFLAGS} -c ${INCLUDE} $<

$(EXE):	${OBJS_MAKE}
	${CC}  ${CFLAGS} -o $@ ${OBJS_MAKE} ${LDFLAGS} -DPAR

debug: ${DOBJS}
	$(CC) -g  ${DFLAGS} -o $@ ${DOBJS} $(LIBS) -DPAR

clean_debug: 
	rm -f $(DEBUGFOLDER)/*.o
	rm -f debug

clean:
	rm -f *.o a.out debug.exe
	rm -f $(MAINOBJFOLDER)/*.o
clear: 
	rm -f *.bin *.dat debug.exe *.all


