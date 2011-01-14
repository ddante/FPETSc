ALL:ex1f

# PETSC_FC_INCLUDES MUST be used in the make file 
# if you don't want to use the PETSc's rules for compiling
# (i.e. include ${PETSC_DIR}/conf/rules)
FFLAGS =  ${PETSC_FC_INCLUDES} -fpp

# Only variables are included, not the rules
include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules

# Linking the the local and the PETSc's libraries
ex1f:ex1f.o
	ifort -o ex1f ex1f.o ${PETSC_LIB}
