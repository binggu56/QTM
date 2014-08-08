QTM
===

quantum trajectory code using mpi 

################################################
QTM code with friction MPI version
###############################################

Description:

      Rcutoff of short-range and long-range is 2.2*Rnn.

History:

Ganrantee corrent version : 1.3.2, 1.4.0 (without long range)


## 1.4.8 

- take advantage of modules of Fortran90 
- put things which is not chaning after initilization into modules 


## 1.4.7 

- optimize the code and put some subroutines into module
- add subroutines to compute correction function 
- instead of broadcast input variables, every processor reads input file 


############################################
	1.4.4 
###################################

1. change mass and density to see the effects on PDF

2. update long-range correction 



#################################
	1.4.3

1. version on darter 




##################################
	1.4.1

1. add long-range correction

2. record trajectories at the end


##########################
      1.4.0

1. parallel the part for computing approximated p,r
   whose size is of ndim*ndim*ntraj

2. fix a bug (didn't get valuse for cp2,cr2) that exist in lower version than 1.3.2

###################################
      1.3.2 

Date: 5/5/2014

Changes:

      1. exclude the long-term potential

####################################
      1.3.1
####################################
1. make Ntraj as a variable in INPUT

2. deep parallel in fitting procedure

3. great improvement in computing time
######################################
      1.2.0
####################################
1. parallel quantum potential part using MPI  

#######################
      1.1.0
######################

1. instead of MPI_BCAST & MPI_REDUCE to send message, use 
   MPI_SCATTER & MPI_GATHER to reduce communication  time

###############################
      1.0.2
##############################
1. add checkpoints  

2. add option to read input or not




