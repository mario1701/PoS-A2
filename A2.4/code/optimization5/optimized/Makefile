####################################################################
#                                                                  #
#             Makefile for FIRE Solver Benchmarks                  #
#                                                                  #
####################################################################

CC = $(PREP) mpicc
CFLAGS = -Wall -g -O3 -mavx $(METIS_INC) $(PAPI_INC) # -lmetis 
LIBS = -lm $(METIS_LIB) $(PAPI_LIB)

LIBPOS=libpos.a
AR = ar
ARFLAGS = rv

SRCS = initialization.c compute_solution.c finalization.c util_read_files.c util_write_files.c domain_distribution.c
OBJS =  $(addsuffix .o, $(basename $(SRCS)))

all: gccg 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

gccg: gccg.c test_functions.c $(LIBPOS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(LIBPOS) : $(OBJS)
	$(AR) $(ARFLAGS) $@ $+

clean:
	rm -rf *.o gccg $(LIBPOS)
