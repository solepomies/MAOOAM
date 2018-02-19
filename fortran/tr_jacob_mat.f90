
! tr_jacob_mat.f90
!
!> Tests to obtain the trace of the Jacobian matrix
!     
!> @copyright                                                               
!> 2016 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!


PROGRAM tr_jacob_mat
  USE params, only:ndim
  USE aotensor_def, only: init_aotensor
  USE tl_ad_tensor, only: init_tltensor, init_adtensor, jacobian_mat
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: y0
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: J
  REAL(KIND=8) :: tr
  INTEGER :: i
  
  ! Compute the tensors

  CALL init_aotensor
  CALL init_tltensor
  CALL init_adtensor


  ALLOCATE(y0(0:ndim),J(ndim,ndim))
  y0=1.D0
  J=jacobian_mat(y0)
  tr=0.D0
  
  DO i=1,ndim
    tr=tr+J(i,i)
  ENDDO

  print*, 'Tr(J) = ',tr
 


END PROGRAM tr_jacob_mat

