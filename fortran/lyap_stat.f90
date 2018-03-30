
! lyap_stat.f90
!
!>  Statistics accumulators for the Lyapunov exponents
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!



MODULE lyap_stat
  USE params, only: ndim
  IMPLICIT NONE

  PRIVATE
  
  INTEGER :: i=0 !< Number of stats accumulated
  
  ! Vectors holding the stats
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m       !< Vector storing the inline mean
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: mprev   !< Previous mean vector
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: v       !< Vector storing the inline variance
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: mtmp  


  PUBLIC :: lyap_acc,lyap_init_stat,lyap_mean,lyap_var,lyap_iter,lyap_reset

  CONTAINS

    !> Initialise the accumulators
    SUBROUTINE lyap_init_stat
      INTEGER :: AllocStat
      
      ALLOCATE(m(0:ndim),mprev(0:ndim),v(0:ndim),mtmp(0:ndim), STAT=AllocStat)
      IF (AllocStat /= 0) STOP '*** Not enough memory ***'
      m=0.D0
      mprev=0.D0
      v=0.D0
      mtmp=0.D0
      
    END SUBROUTINE lyap_init_stat

    !> Accumulate one state
    SUBROUTINE lyap_acc(x)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: x
      i=i+1
      mprev=m+(x-m)/i
      mtmp=mprev
      mprev=m
      m=mtmp
      v=v+(x-mprev)*(x-m)
    END SUBROUTINE lyap_acc

    !> Function returning the mean
    FUNCTION lyap_mean()
      REAL(KIND=8), DIMENSION(0:ndim) :: lyap_mean
      lyap_mean=m
    END FUNCTION lyap_mean

    !> Function returning the variance
    FUNCTION lyap_var()
      REAL(KIND=8), DIMENSION(0:ndim) :: lyap_var
      lyap_var=v/(i-1)
    END FUNCTION lyap_var

    !> Function returning the number of data accumulated
    FUNCTION lyap_iter()
      INTEGER :: lyap_iter
      lyap_iter=i
    END FUNCTION lyap_iter

    !> Routine resetting the accumulators
    SUBROUTINE lyap_reset
      m=0.D0
      mprev=0.D0
      v=0.D0
      i=0
    END SUBROUTINE lyap_reset
      

  END MODULE lyap_stat
