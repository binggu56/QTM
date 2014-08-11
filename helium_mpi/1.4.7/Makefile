objects = mod.o long_range.o reduce.o qm.o derivs.o input.o tstamp.o vinit.o gasdev.o ran1.o fit.o corr.o mativs.o
#objects2 = long_range.o reduce.o cont.o derivs.o input.o tstamp.o vinit.o qpot.o gasdev.o ran2.o
FC = mpif90
#FC = ftn
switch = -O3 #-w
EXEC = qm.x
#FLAGS = -heap-arrays
#EXEC2 = cont
#LIB = /usr/lib64/atlas/liblapack.so.3.0
LIB = -llapack  
#LIB = -mkl=parallel

$(EXEC): $(objects)
	$(FC) -o $(EXEC) $(FLAGS) $(LIB) $(switch) $(objects)
#$(EXEC2): $(objects2)
#	$(FC) -o $(EXEC2) $(LIB) $(switch) $(objects2)

%.mod: %.f90
	$(FC) -c $<

%.o: %.f
	$(FC) -c $<

%.o: %.f90
	$(FC) -c $<

clean:
	rm $(objects) $(EXEC)
veryclean:
	rm *.o *.dat *.mod
