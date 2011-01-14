ALL: test

FFLAGS =  ${PETSC_FC_INCLUDES} -fpp
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

test: solv_lin_syst.o main.o 
	ifort -o test solv_lin_syst.o main.o ${PETSC_LIB}


