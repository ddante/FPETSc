MODULE solv_lin_syst

IMPLICIT NONE

CONTAINS

   !-------------------------------------
   SUBROUTINE PETSc_lin_solv(V, I, J, nf)
   !-------------------------------------
   !
   IMPLICIT NONE

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"

!  Variables:
!     ksp     - linear solver context
!     ksp      - Krylov subspace method context
!     pc       - preconditioner context
!     x, b, u  - approx solution, right-hand-side, exact solution vectors
!     A        - matrix that defines linear system
!     its      - iterations for convergence
!     norm     - norm of error in solution
!
      Vec              x, b, u
      Mat              A 
      KSP              ksp
      PC               pc
      PetscReal        norm, tol
      PetscErrorCode   ierr
      PetscTruth       flg
      PetscMPIInt      size,rank     
      
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: V
      INTEGER,      DIMENSION(:), INTENT(IN) :: I
      INTEGER,      DIMENSION(:), INTENT(IN) :: J
      INTEGER,                    INTENT(IN) :: nf
      !=============================================

      PetscInt,DIMENSION(nf)  :: nz
      PetscInt  :: n, n1, row0, its, max_res
      PetscReal :: one, none
      PetscInt,  DIMENSION(:), ALLOCATABLE :: col
      PetscReal, DIMENSION(:), ALLOCATABLE :: Mval
      
      INTEGER :: row_id, k      
      
      n = nf
      n1 = 1
      one = 1.0
      none = -1.0
      
      !-------------------------
      ! INIZIALIZATION of PETSc
      !------------------------------------------------      
      CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
      
      CALL MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)
      IF (size /= 1) THEN
         CALL MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
         IF (rank == 0) THEN
            WRITE(6,*) 'This is a uniprocessor example only!'
         ENDIF
         SETERRQ(1, ' ', ierr)
      ENDIF

      CALL PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-n', n, flg, ierr)
      
      write(*,*) '*** PETSc initializated ***'
      
      !----------------------------------------------
      ! CREATION & ASSEMBLE of the MATRIX and VECTORS
      !------------------------------------------------------- 
      
      DO k = 1, n
         nz(k) = I(k+1) - I(k)
      ENDDO
        
      CALL MatCreateSeqAIJ(PETSC_COMM_SELF, n,n, PETSC_NULL_INTEGER, nz, A, ierr)
      
      WRITE(*,*) '*** Matrix A initializated ***'
           
      DO k = 1, n
      
          ALLOCATE(col(nz(k)), Mval(nz(k)))
              
         row_id = I(k)
           Mval = V( row_id : (row_id+nz(k)-1) )
            col = J( row_id : (row_id+nz(k)-1) )
         
         row0 =   k-1
          col = col-1
           
         CALL MatSetValues(A, n1,row0, nz(k),col, Mval, INSERT_VALUES, ierr)
      
         DEALLOCATE(col, Mval)
      
      ENDDO
      
      CALL MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
      CALL MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY, ierr)

      WRITE(*,*) '*** Matrix A assembled ***'
      
      !!CALL MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr)      
      
      ! Create solution vector
      CALL VecCreate(PETSC_COMM_WORLD, x, ierr)
      CALL VecSetSizes(x, PETSC_DECIDE, n, ierr)
      CALL VecSetFromOptions(x, ierr)
      CALL VecDuplicate(x, u, ierr)
      
      CALL VecSet(u, one, ierr)
      
      WRITE(*,*) '*** Vectors of the solution created ***'
      
      ! Create RHS vector
      CALL VecCreate(PETSC_COMM_WORLD, b, ierr)
      CALL VecSetSizes(b, PETSC_DECIDE, n, ierr)
      CALL VecSetFromOptions(b, ierr)
      
      CALL MatMult(A,u, b, ierr)
      
      WRITE(*,*) '*** RHS vector created ***'
      
      !----------------------------------------------
      ! SOLVE the LINEAR SYSTEM
      !------------------------------------------------------- 
      CALL KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
      
      WRITE(*,*) '*** KSP context created ***' 

      CALL KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN,ierr)
      
      WRITE(*,*) '*** KSP operators set ***'  
      
      CALL KSPSetType(ksp, KSPGMRES, ierr)
      max_res = 30
      CALL KSPGMRESSetRestart(ksp, max_res, ierr)
      
      ! PRECONDITIONING
      CALL KSPGetPC(ksp, pc, ierr)
      CALL PCSetType(pc, PCLU, ierr)
      
      tol = 1.d-14
      
      CALL KSPSetTolerances(ksp, tol, PETSC_DEFAULT_DOUBLE_PRECISION,     &
     &     PETSC_DEFAULT_DOUBLE_PRECISION, 100000, ierr)
            
      CALL KSPSetFromOptions(ksp, ierr)
      
      CALL KSPSetInitialGuessNonzero(ksp, flg, ierr)
      CALL KSPSolve(ksp, b, x, ierr)

      !CALL VecView(x, PETSC_VIEWER_STDOUT_WORLD, ierr)

      CALL KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD, ierr)  
      
      WRITE(*,*) '*** LINEAR SYSTEM SOLVED ***'
      
     
! !!!!!!!!!!!!!!!!!!!     
      CALL VecAXPY(x, none, u, ierr)
      CALL VecNorm(x, NORM_2, norm, ierr)
      CALL KSPGetIterationNumber(ksp, its, ierr)
      IF (norm >= 1.e-12) THEN
        WRITE(6,100) norm,its
      ELSE
        WRITE(6,200) its
      ENDIF
      
 100  FORMAT('Norm of error = ',e12.4,',  Iterations = ',i5)
 200  FORMAT('Norm of error < 1.e-12,Iterations = ',i5)
! !!!!!!!!!!!!!!!!!!!

      !----------
      ! CLEAN UP
      !--------------------------- 

      CALL VecDestroy(x,ierr)
      CALL VecDestroy(u,ierr)
      CALL VecDestroy(b,ierr)
      CALL MatDestroy(A,ierr)
      CALL KSPDestroy(ksp,ierr)
      CALL PetscFinalize(ierr)
            
    END SUBROUTINE PETSc_lin_solv
   !----------------------------

END MODULE solv_lin_syst
