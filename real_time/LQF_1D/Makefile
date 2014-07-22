objects = qm.o ran2.o mod.o qpot.o fcl.o gasdev.o
FC = gfortran
switch = -O3
EXEC = qm
LIB = /usr/lib64/atlas/liblapack.so.3.0
#LIB = -llapack  
$(EXEC): $(objects)
	$(FC) -o $(EXEC) $(LIB) $(switch) $(objects)

cdat.mod:mod.o mod.f
	$(FC) -c mod.f
gauss.o:cdat.mod gauss.f
	$(FC) -c gauss.f

%.o: %.f 
	$(FC) -c  $<
#%.O: %.f
#	$(FC) -c $<

clean:
	rm $(objects)
veryclean:
	rm *.o
