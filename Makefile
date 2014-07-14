objects = mod.o long_range.o reduce.o qm.o derivs.o input.o tstamp.o vinit.o gasdev.o ran1.o fit.o
#objects2 = long_range.o reduce.o cont.o derivs.o input.o tstamp.o vinit.o qpot.o gasdev.o ran2.o
#FC = gfortran
FC = ftn
switch = -O3 #-w
EXEC = qm
#FLAGS = -heap-arrays
#EXEC2 = cont
#LIB = /usr/lib64/atlas/liblapack.so.3.0
#LIB = -llapack  
#LIB = -mkl=sequential 
LIB = -mkl
$(EXEC): $(objects)
	$(FC) -o $(EXEC) $(FLAGS) $(LIB) $(switch) $(objects)
#$(EXEC2): $(objects2)
#	$(FC) -o $(EXEC2) $(LIB) $(switch) $(objects2)

cdat.mod: mod.f
	$(FC) -c mod.f

%.o: %.f
	$(FC) -c $<

%.o: %.f90
	$(FC) -c $<

clean:
	rm $(objects) $(EXEC)
veryclean:
	rm *.o *.dat *.mod
