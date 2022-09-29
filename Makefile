F95 = mpif90
OPTS = -O3 -Wall 
LIBS=-L${BLASDIR} ${BLASLIB} -lpthread

OBJS = header.o initial.o matmult.o rhs.o timestepping.o grayscott.o save_fields.o

all: grayscott

header.mod: header.f90
	$(F95) $(OPT) -c header.f90

%.o:  %.f90 header.mod
	$(F95) $(OPTS) -c $<

grayscott: $(OBJS)
	$(F95) $(OPTS) $(OBJS) -o grayscott $(LIBS)

clean:
	rm -rf grayscott *.o *.mod core.* solution.* slurm-* *~

