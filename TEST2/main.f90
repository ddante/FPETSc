PROGRAM main

USE solv_lin_syst

IMPLICIT NONE

   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: V
   INTEGER, DIMENSION(:), ALLOCATABLE :: I
   INTEGER, DIMENSION(:), ALLOCATABLE :: J
   INTEGER :: n
   INTEGER :: nnze
   
   INTEGER :: k

   n = 100
   
   nnze = (n-2)*3 + 4
   
   ALLOCATE(V(nnze), J(nnze), I(n+1) )
   
   V(1:2) = (/ 2.0, -1.0 /)
   DO k = 1, n-2
       V(k*3: 2+k*3) = (/ -1.0, 2.0, -1.0 /)
   ENDDO      
   V(nnze-1:nnze) = (/ -1.0, 2.0 /)

   
   J(1:2) = (/ 1, 2 /)
   DO k = 1, n-2
       J(k*3: 2+k*3) = (/ k, k+1, k+2 /)
   ENDDO
   J(nnze-1:nnze) =  (/ n-1, n /)
   
   I(1) = 1
   DO k = 1, n-1
       I(k+1) = k*3
   ENDDO
   I(n+1) = nnze+1
   
   CALL PETSc_lin_solv(V, I, J, n)
   

END PROGRAM main
