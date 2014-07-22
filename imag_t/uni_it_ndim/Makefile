objects = qm.o fcl.o qpot.o gasdev.o ran2.o  
objects2 = qm.o fcl.o qpot_linear.o  gasdev.o ran2.o

FC = gfortran
switch = -O3 
EXEC = qm
EXEC2 = qm2
LIB = /usr/lib64/atlas/liblapack.so.3.0
#LIB = -llapack  
FLAGS = -w

$(EXEC): $(objects)
	$(FC) -o $(EXEC) $(LIB) $(FLAGS) $(switch) $(objects)
$(EXEC2): $(objects2)
	$(FC) -o $(EXEC2) $(LIB) $(FLAGS) $(switch) $(objects2)

%.o: %.f
	$(FC) -c  $<
#%.O: %.f
#	$(FC) -c $<

clean:
	rm $(objects)
veryclean:
	rm *.o
