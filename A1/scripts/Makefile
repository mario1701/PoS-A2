####################################################################
#                                                                  #
#             Makefile for FIRE Solver Benchmarks                  #
#                                                                  #
####################################################################

#Compilers
CC = mpicc -icc
CFLAGS = -Wall -g -O0 -vec

#linker
LIBS = -lm $(PAPI_INC) $(PAPI_LIB)

LIBPOS=libpos.a
AR = ar
ARFLAGS = rv

#Include Files
SRCS = initialization.c compute_solution.c finalization.c util_read_files.c util_write_files.c vol2mesh.c
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

all: gccg 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

gccg: gccg.c $(LIBPOS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(LIBPOS) : $(OBJS)
	$(AR) $(ARFLAGS) $@ $+

binconv:
	$(CC) $(CFLAGS) binconv.c -o binconv

clean:
	rm -rf *.o gccg binconv $(LIBPOS)
	

