
! aotensor_def.f90
!
!>  The equation tensor for the coupled ocean-atmosphere model
!>  with temperature which allows for an extensible set of modes
!>  in the ocean and in the atmosphere.
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!> @remark                                                                 
!> Generated Fortran90/95 code 
!> from aotensor.lua
!                                                                           
!---------------------------------------------------------------------------!


MODULE aotensor_def

  !-----------------------------------------------------!
  !                                                     !
  ! Preamble and variables declaration                  !
  !                                                     !
  !-----------------------------------------------------!

  USE params
  USE inprod_analytic
  USE tensor, only:coolist,simplify,add_to_tensor,add_check
  IMPLICIT NONE

  PRIVATE

  !> Vector used to count the tensor elements
  INTEGER, DIMENSION(:), ALLOCATABLE :: count_elems

  !> Epsilon to test equality with 0
  REAL(KIND=8), PARAMETER :: real_eps = 2.2204460492503131e-16

  PUBLIC :: init_aotensor

  !> \f$\mathcal{T}_{i,j,k}\f$ - Tensor representation of the tendencies.
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: aotensor

  !> Buffer for the aotensor calculation
  TYPE(coolist), DIMENSION(:), ALLOCATABLE :: aobuf


  !-----------------------------------------------------!
  !                                                     !
  ! End of preamble                                     !
  !                                                     !
  !-----------------------------------------------------!

CONTAINS

  !-----------------------------------------------------!
  !                                                     !
  ! Function declarations                               !
  !                                                     !
  !-----------------------------------------------------!

  !> Translate the \f$\psi_{a,i}\f$ coefficients into effective coordinates
  FUNCTION psi(i)
    INTEGER :: i,psi
    psi = i
  END FUNCTION psi

  !> Translate the \f$\theta_{a,i}\f$ coefficients into effective coordinates
  FUNCTION theta(i)
    INTEGER :: i,theta
    theta = i + natm
  END FUNCTION theta

  !> Translate the \f$\psi_{o,i}\f$ coefficients into effective coordinates
  FUNCTION A(i)
    INTEGER :: i,A
    A = i + 2 * natm
  END FUNCTION A

  !> Translate the \f$\delta T_{o,i}\f$ coefficients into effective coordinates
  FUNCTION T(i)
    INTEGER :: i,T
    T = i + 2 * natm + noc
  END FUNCTION T

  !> Kronecker delta function
  FUNCTION kdelta(i,j)
    INTEGER :: i,j,kdelta
    kdelta=0
    IF (i == j) kdelta = 1
  END FUNCTION kdelta


  !> Subroutine to compute the tensor #aotensor
  SUBROUTINE compute_aotensor
    INTEGER :: i,j,k,l,n
    REAL(KIND=8) :: v
    CALL add_check(aobuf,theta(1),0,0,(Cpa / (1 - atmos%a(1,1) * sig0)),aotensor)
    DO i = 1, natm
       DO j = 1, natm
          CALL add_check(aobuf,psi(i),psi(j),0,-(((atmos%c(i,j) * betp) / atmos%a(i,i))) - (kd * kdelta(i,j)) / 2,aotensor)
          CALL add_check(aobuf,theta(i),psi(j),0,(atmos%a(i,j) * kd * sig0) / (-2 + 2 * atmos%a(i,i) * sig0),aotensor)
          CALL add_check(aobuf,psi(i),theta(j),0,(kd * kdelta(i,j)) / 2,aotensor)
          CALL add_check(aobuf,theta(i),theta(j),0,(-((sig0 * (2. * atmos%c(i,j) * betp +&
               & atmos%a(i,j) * (kd + 4. * kdp)))) + 2. * (LSBpa + sc * Lpa) &
               &* kdelta(i,j)) / (-2. + 2. * atmos%a(i,i) * sig0),aotensor)
       END DO
    END DO
    DO i = 1, natm
       DO n=1,atmos%g(i)%nelems
          j=atmos%g(i)%elems(n)%j
          k=atmos%g(i)%elems(n)%k
          v=atmos%g(i)%elems(n)%v
          CALL add_check(aobuf,psi(i),psi(j),psi(k),-((v*atmos%a(k,k) / atmos%a(i,i))),aotensor)
          CALL add_check(aobuf,psi(i),theta(j),theta(k),-((v*atmos%a(k,k) / atmos%a(i,i))),aotensor)
          CALL add_check(aobuf,theta(i),psi(j),theta(k),(v*(1- atmos%a(k,k)* sig0)  / (-1 + atmos%a(i,i)*sig0)),aotensor)
          CALL add_check(aobuf,theta(i),theta(j),psi(k),(v*atmos%a(k,k) * sig0) / (1 - atmos%a(i,i) * sig0),aotensor)
       ENDDO
    END DO
    DO i = 1, natm
       DO j = 1, noc
          CALL add_check(aobuf,psi(i),A(j),0,kd * atmos%d(i,j) / (2 * atmos%a(i,i)),aotensor)
          CALL add_check(aobuf,theta(i),A(j),0,kd * (atmos%d(i,j) * sig0) / (2 - 2 * atmos%a(i,i) * sig0),aotensor)
          CALL add_check(aobuf,theta(i),T(j),0,atmos%s(i,j) * (2 * LSBpo + Lpa) / (2 - 2 * atmos%a(i,i) * sig0),aotensor)
       END DO
    END DO
    DO i = 1, noc
       DO j = 1, natm
          CALL add_check(aobuf,A(i),psi(j),0,ocean%K(i,j) * dp / (ocean%M(i,i) + G),aotensor)
          CALL add_check(aobuf,A(i),theta(j),0,-(ocean%K(i,j)) * dp / (ocean%M(i,i) + G),aotensor)
       END DO
    END DO
    DO i = 1, noc
       DO j = 1, noc
          CALL add_check(aobuf,A(i),A(j),0,-((ocean%N(i,j) * betp + ocean%M(i&
               &,i) * (rp + dp) * kdelta(i,j))) / (ocean%M(i,i) + G),aotensor)
       END DO
    END DO
    DO i = 1, noc
       DO n=1,ocean%O(i)%nelems
          j=ocean%O(i)%elems(n)%j
          k=ocean%O(i)%elems(n)%k
          v=ocean%O(i)%elems(n)%v
          CALL add_check(aobuf,A(i),A(j),A(k),-(v*ocean%M(k,k)) / (ocean%M(i,i) + G),aotensor)
       END DO
    END DO
    DO i = 1, noc
       CALL add_check(aobuf,T(i),0,0,Cpo * ocean%W(i,1),aotensor)
       DO j = 1, natm
          CALL add_check(aobuf,T(i),theta(j),0,ocean%W(i,j) * (2 * sc * Lpo + sBpa),aotensor)
       END DO
    END DO
    DO i = 1, noc
       DO j = 1, noc
          CALL add_check(aobuf,T(i),T(j),0,-((Lpo + sBpo)) * kdelta(i,j),aotensor)
       END DO
    END DO
    DO i = 1, noc
       DO n=1,ocean%O(i)%nelems
          j=ocean%O(i)%elems(n)%j
          k=ocean%O(i)%elems(n)%k
          v=ocean%O(i)%elems(n)%v
          CALL add_check(aobuf,T(i),A(j),T(k),-v,aotensor)
       END DO
    END DO
    CALL add_to_tensor(aobuf,aotensor)
  END SUBROUTINE compute_aotensor

  !-----------------------------------------------------!
  !                                                     !
  ! Initialisation routine                              !
  !                                                     !
  !-----------------------------------------------------!

  !> Subroutine to initialise the #aotensor tensor
  !> @remark This procedure will also call params::init_params() and inprod_analytic::init_inprod() .
  SUBROUTINE init_aotensor
    INTEGER :: i
    INTEGER :: AllocStat 

    CALL init_params  ! Iniatialise the parameter

    CALL init_inprod  ! Initialise the inner product tensors

    ALLOCATE(aotensor(ndim),aobuf(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    DO i=1,ndim
       ALLOCATE(aobuf(i)%elems(1000), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    END DO
    
    CALL compute_aotensor

    DEALLOCATE(aobuf, STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"

    CALL simplify(aotensor)

  END SUBROUTINE init_aotensor
END MODULE aotensor_def
      


