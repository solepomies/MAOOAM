
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

!  USE params
  USE params
  USE inprod_analytic
  USE tensor, only:coolist,simplify
  IMPLICIT NONE

  PRIVATE

  !> Vector used to count the tensor elements
  INTEGER, DIMENSION(:), ALLOCATABLE :: count_elems

  !> Epsilon to test equality with 0
  REAL(KIND=8), PARAMETER :: real_eps = 2.2204460492503131e-16

  PUBLIC :: init_aotensor

  !> \f$\mathcal{T}_{i,j,k}\f$ - Tensor representation of the tendencies.
  TYPE(coolist), DIMENSION(:), ALLOCATABLE, PUBLIC :: aotensor


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

  !> Subroutine to add element in the #aotensor \f$\mathcal{T}_{i,j,k}\f$ structure.
  !> @param i tensor \f$i\f$ index
  !> @param j tensor \f$j\f$ index
  !> @param k tensor \f$k\f$ index
  !> @param v value to add
  SUBROUTINE coeff(i,j,k,v)
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN) :: v
    INTEGER :: n
    IF (.NOT. ALLOCATED(aotensor)) STOP "*** coeff routine : tensor not yet allocated ***"
    IF (.NOT. ALLOCATED(aotensor(i)%elems)) STOP "*** coeff routine : tensor not yet allocated ***"
    IF (ABS(v) .ge. real_eps) THEN
       n=(aotensor(i)%nelems)+1
       IF (j .LE. k) THEN
          aotensor(i)%elems(n)%j=j
          aotensor(i)%elems(n)%k=k
       ELSE
          aotensor(i)%elems(n)%j=k
          aotensor(i)%elems(n)%k=j
       END IF
       aotensor(i)%elems(n)%v=v
       aotensor(i)%nelems=n
    END IF
  END SUBROUTINE coeff

  !> Subroutine to count the elements of the #aotensor \f$\mathcal{T}_{i,j,k}\f$. Add +1 to count_elems(i) for each value that is added to the tensor i-th component.
  !> @param i tensor \f$i\f$ index
  !> @param j tensor \f$j\f$ index
  !> @param k tensor \f$k\f$ index
  !> @param v value that will be added
  SUBROUTINE add_count(i,j,k,v)
    INTEGER, INTENT(IN) :: i,j,k
    REAL(KIND=8), INTENT(IN)  :: v
    IF (ABS(v) .ge. real_eps) count_elems(i)=count_elems(i)+1
  END SUBROUTINE add_count

  !> Subroutine to compute the tensor #aotensor
  !> @param func External function to be used
  SUBROUTINE compute_aotensor(func)
    EXTERNAL :: func
    INTERFACE
       SUBROUTINE func(i,j,k,v)
         INTEGER, INTENT(IN) :: i,j,k
         REAL(KIND=8), INTENT(IN) :: v
       END SUBROUTINE func
    END INTERFACE
    INTEGER :: i,j,k
    CALL func(theta(1),0,0,(Cpa / (1 - atmos%a(1,1) * sig0))) ! ok
    CALL func(T(1),0,0,Cpo) ! je l'ai déplacé car il n'est plus non nul qu'en i = 1 ; les autres termes de la même "ligne" seront-ils bien initialisés à 0, ou a-t-on vraiment besoin d'une itération sur j avec la fonction de Kronecker ?
    DO i = 1, natm
       DO j = 1, natm
          CALL func(psi(i),psi(j),0,-(((atmos%c(correspatm(i),correspatm(j)) * betp) / atmos%a(correspatm(i),correspatm(i)))) -&
               &(kd * kdelta(i,j)) / 2 + atmos%a(i,j)*nuap) ! le dernier terme est un terme de dissipation qui n'est pas mentionné dans le papier
          CALL func(theta(i),psi(j),0,(atmos%a(correspatm(i),correspatm(j)) * kd * sig0) / (-2 + 2 * atmos%a(correspatm(i),&
               &correspatm(i)) * sig0)) ! ok
          CALL func(psi(i),theta(j),0,(kd * kdelta(i,j)) / 2) ! ok
          CALL func(theta(i),theta(j),0,(-(sig0 * (2. * atmos%c(correspatm(i),correspatm(j)) * betp +&
               & atmos%a(correspatm(i),correspatm(j)) * (kd + 4. * kdp))) + 2. * (LSBpa + sc * Lpa) &
               &* kdelta(i,j)) / (-2. + 2. * atmos%a(correspatm(i),correspatm(i)) * sig0)) ! le facteur sc vient de la différence entre la température moyenne de l'atmosphère et sa température à la surface de l'océan : dans l'équation d'évolution de la température de l'atmosphère, on écrira plutôt le terme proportionnel à la différence des températures comme étant proportionnel à la différence entre la température océanique et sc*(la température atmosphérique)
          DO k = 1, natm
             CALL func(psi(i),psi(j),psi(k),-((atmos%b(correspatm(i),correspatm(j),correspatm(k)) / atmos%a(correspatm(i),&
                  &correspatm(i))))) ! ok
             CALL func(psi(i),theta(j),theta(k),-((atmos%b(correspatm(i),correspatm(j),correspatm(k)) / atmos%a(correspatm(i),&
                  &correspatm(i))))) ! ok
             CALL func(theta(i),psi(j),theta(k),(atmos%g(correspatm(i),correspatm(j),correspatm(k)) -&
                  & atmos%b(correspatm(i),correspatm(j),correspatm(k)) * sig0) / (-1 + atmos%a(correspatm(i),correspatm(i)) *&
                  & sig0)) ! ok
             CALL func(theta(i),theta(j),psi(k),(atmos%b(correspatm(i),correspatm(j),correspatm(k)) * sig0) / (1 -&
                  & atmos%a(correspatm(i),correspatm(i)) * sig0)) ! ok
          END DO
       END DO
       DO j = 1, noc
          CALL func(psi(i),A(j),0,kd * atmos%a(correspatm(i),correspoc(j)) / (2 * atmos%a(correspatm(i),correspatm(i)))) ! ok
          CALL func(theta(i),A(j),0,kd * (atmos%a(correspatm(i),correspoc(j)) * sig0) / (2 - 2 * atmos%a(correspatm(i),&
               &correspatm(i)) * sig0)) ! ok
          CALL func(theta(i),T(j),0,kdelta(i,j) * (2 * LSBpo + Lpa) / (2 - 2 * atmos%a(correspatm(i),correspatm(i)) * sig0)) ! ok
       END DO
    END DO
    DO i = 1, noc
       DO j = 1, natm
          CALL func(A(i),psi(j),0,atmos%a(correspoc(i),correspatm(j)) * dp / (atmos%a(correspoc(i),correspoc(i)) + G)) ! ok
          CALL func(A(i),theta(j),0,-(atmos%a(correspoc(i),correspatm(j))) * dp / (atmos%a(correspoc(i),correspoc(i)) + G)) ! ok
       END DO
       DO j = 1, noc
          CALL func(A(i),A(j),0,-((atmos%c(correspoc(i),correspoc(j)) * betp + atmos%a(correspoc(i),correspoc(i)) * (rp + dp) *&
               & kdelta(i,j) - atmos%a(correspoc(i),correspoc(j))**2*nuop)) / (atmos%a(correspoc(i),correspoc(i)) + G)) ! le dernier terme est un terme de dissipation qui n'est pas mentionné dans le papier
          DO k = 1, noc
             CALL func(A(i),A(j),A(k),-(atmos%b(correspoc(i),correspoc(j),correspoc(k))) / (atmos%a(correspoc(i),correspoc(i)) +&
                  & G)) ! ok
          END DO
       END DO
    END DO
    DO i = 1, noc       
       CALL func(T(i),theta(i),0,(2 * sc * Lpo + sBpa)) ! le facteur sc vient de la modification de l'équation d'évolution de la température atmosphérique mentionnée plus haut
       DO j = 1, noc
          CALL func(T(i),T(j),0,-((Lpo + sBpo)) * kdelta(i,j)) ! ok
          DO k = 1, noc
             CALL func(T(i),A(j),T(k),-(atmos%g(correspoc(i),correspoc(j),correspoc(k)))) ! ok
          END DO
       END DO
    END DO
  END SUBROUTINE compute_aotensor

  !-----------------------------------------------------!
  !                                                     !
  ! Initialisation routine                              !
  !                                                     !
  !-----------------------------------------------------!

  !> Subroutine to initialise the #aotensor tensor
  !> @remark This procedure will also call params::init_params() and inprod_analytic::init_inprod() .
  !> It will finally call inprod_analytic::deallocate_inprod() to remove the inner products, which are not needed
  !> anymore at this point.
  SUBROUTINE init_aotensor
    INTEGER :: i
    INTEGER :: AllocStat 

    CALL init_params  ! Iniatialise the parameter

    CALL init_inprod  ! Initialise the inner product tensors

    ALLOCATE(aotensor(ndim),count_elems(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    count_elems=0

    CALL compute_aotensor(add_count)

    DO i=1,ndim
       ALLOCATE(aotensor(i)%elems(count_elems(i)), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    END DO

    DEALLOCATE(count_elems, STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"
    
    CALL compute_aotensor(coeff)

    CALL simplify(aotensor)

  END SUBROUTINE init_aotensor
END MODULE aotensor_def
      


