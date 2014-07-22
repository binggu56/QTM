# Simple makefile
.SUFFIXES: .f .x
.f.x:
	$(FC) $(Flags) -ffixed-form $(LIB) $< -o $*.x $(SRC)
  
    FC = gfortran
#    FC = /opt/intel/composer_xe_2011_sp1.8.273/bin/ia32/ifort

#  Flags =  -O0 -ffixed-line-length-none 
#    Flags = -llapack 

#  LIB = /opt/intel/composer_xe_2011_sp1.8.273/mkl/interfaces/lapack95
   LIB = /usr/lib64/atlas/liblapack.so.3
   SRC = random.f gasdev.f ran2.f quantum_potential.f fcl.f
