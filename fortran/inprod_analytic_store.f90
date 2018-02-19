
! inprod_analytic.f90
!
!>  Inner products between the truncated set of basis functions for the
!>  ocean and atmosphere streamfunction fields.
!>  These are partly calculated using the analytical expressions from
!>  Cehelsky, P., & Tung, K. K. : Theories of multiple equilibria and
!>  weather regimes-A critical reexamination. Part II: Baroclinic two-layer
!>  models. Journal of the atmospheric sciences, 44(21), 3282-3303, 1987.
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!> @remark                                                                 
!> Generated Fortran90/95 code 
!> from inprod_analytic.lua
!                                                                           
!---------------------------------------------------------------------------!

MODULE inprod_analytic

  !-----------------------------------------------------!
  !                                                     !
  ! Preamble and variables declaration                  !
  !                                                     !
  !-----------------------------------------------------!

  USE params, only: nbatm, nboc, natm, noc, ndim, n, oms, ams, pi
  USE tensor, only:coolist,add_check,add_to_tensor,load_tensor_from_file,write_tensor_to_file
  USE util, only: isin,piksrt
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_inprod

  !> Atmospheric bloc specification type
  TYPE :: atm_wavenum 
     CHARACTER :: typ
     INTEGER :: M=0,P=0,H=0
     REAL(KIND=8) :: Nx=0.,Ny=0.
  END TYPE atm_wavenum

  !> Oceanic bloc specification type
  TYPE :: ocean_wavenum
     INTEGER :: P,H
     REAL(KIND=8) :: Nx,Ny
  END TYPE ocean_wavenum

  !> Type holding the atmospheric inner products tensors
  TYPE :: atm_tensors
     REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: a,c,d,s
     TYPE(coolist), DIMENSION(:), ALLOCATABLE :: g
  END TYPE atm_tensors

  !> Type holding the oceanic inner products tensors
  TYPE :: ocean_tensors
     REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: K,M,N,W
     TYPE(coolist), DIMENSION(:), ALLOCATABLE :: O
  END TYPE ocean_tensors

  !> Atmospheric blocs specification
  TYPE(atm_wavenum), DIMENSION(:), ALLOCATABLE, PUBLIC :: awavenum 
  !> Oceanic blocs specification
  TYPE(ocean_wavenum), DIMENSION(:), ALLOCATABLE, PUBLIC :: owavenum 

  !> Atmospheric tensors
  TYPE(atm_tensors), PUBLIC :: atmos 
  !> Oceanic tensors
  TYPE(ocean_tensors), PUBLIC :: ocean

  !> Buffer for the inner products calculation
  TYPE(coolist), DIMENSION(:), ALLOCATABLE :: ipbuf


  !-----------------------------------------------------!
  !                                                     !
  ! End of preamble                                     !
  !                                                     !
  !-----------------------------------------------------!

CONTAINS

  !-----------------------------------------------------!
  !                                                     !
  ! Definition of the Helper functions from Cehelsky    !
  ! & Tung                                              !
  !                                                     !
  !-----------------------------------------------------!

  !>  Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION B1(Pi, Pj, Pk)
    INTEGER :: Pi,Pj,Pk
    B1 = (Pk + Pj) / REAL(Pi)
  END FUNCTION B1

  !>  Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION B2(Pi, Pj, Pk)
    INTEGER :: Pi,Pj,Pk
    B2 = (Pk - Pj) / REAL(Pi)
  END FUNCTION B2

  !>  Integer Dirac delta function
  REAL(KIND=8) FUNCTION delta(r)
    INTEGER :: r
    IF (r==0) THEN
       delta = 1.D0
    ELSE
       delta = 0.D0
    ENDIF
  END FUNCTION delta

  !>  "Odd or even" function
  REAL(KIND=8) FUNCTION flambda(r)
    INTEGER :: r
    IF (mod(r,2)==0) THEN
       flambda = 0.D0
    ELSE
       flambda = 1.D0
    ENDIF
  END FUNCTION flambda

  !>  Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION S1(Pj, Pk, Mj, Hk)
    INTEGER :: Pk,Pj,Mj,Hk
    S1 = -((Pk * Mj + Pj * Hk)) / 2.D0
  END FUNCTION S1

  !>  Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION S2(Pj, Pk, Mj, Hk)
    INTEGER :: Pk,Pj,Mj,Hk
    S2 = (Pk * Mj - Pj * Hk) / 2.D0
  END FUNCTION S2

  !>  Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION S3(Pj, Pk, Hj, Hk)
    INTEGER :: Pj,Pk,Hj,Hk
    S3 = (Pk * Hj + Pj * Hk) / 2.D0
  END FUNCTION S3

  !>  Cehelsky & Tung Helper functions
  REAL(KIND=8) FUNCTION S4(Pj, Pk, Hj, Hk)
    INTEGER :: Pj,Pk,Hj,Hk
    S4 = (Pk * Hj - Pj * Hk) / 2.D0
  END FUNCTION S4

  !-----------------------------------------------------!
  ! Inner products definition routines                  !
  !--------------------------------------------------------!
  ! 1. Inner products in the equations for the atmosphere  !
  !--------------------------------------------------------!
  
  !> Eigenvalues of the Laplacian (atmospheric)
  !> 
  !> \f$ a_{i,j} = (F_i, \nabla^2 F_j)\f$ .
  SUBROUTINE calculate_a
    INTEGER :: i
    TYPE(atm_wavenum) :: Ti
    INTEGER :: AllocStat 
    IF (natm == 0 ) THEN
       STOP "*** Problem with calculate_a : natm==0 ! ***"
    ELSE
       IF (.NOT. ALLOCATED(atmos%a)) THEN
          ALLOCATE(atmos%a(natm,natm), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       END IF
    END IF
    atmos%a=0.D0

    DO i=1,natm
       Ti = awavenum(i)
       atmos%a(i,i) = -(n**2) * Ti%Nx**2 - Ti%Ny**2
    ENDDO
  END SUBROUTINE calculate_a


  !> Beta term for the atmosphere
  !> 
  !> \f$ c_{i,j} = (F_i, \partial_x F_j)\f$ .
  !
  !> @remark
  !> Strict function !! Only accepts KL type.
  !> For any other combination, it will not calculate anything
  SUBROUTINE calculate_c_atm
    INTEGER :: i,j
    TYPE(atm_wavenum) :: Ti, Tj
    REAL(KIND=8) :: val
    INTEGER :: AllocStat 

    IF (natm == 0 ) THEN
       STOP "*** Problem with calculate_c_atm : natm==0 ! ***"
    ELSE
       IF (.NOT. ALLOCATED(atmos%c)) THEN
          ALLOCATE(atmos%c(natm,natm), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       END IF
    END IF
    atmos%c=0.D0

    DO i=1,natm
       DO j=i,natm
          Ti = awavenum(i)
          Tj = awavenum(j)
          val = 0.D0
          IF ((Ti%typ == "K") .AND. (Tj%typ == "L")) THEN 
             val = n * Ti%M * delta(Ti%M - Tj%H) * delta(Ti%P - Tj%P)
          END IF
          IF (val /= 0.D0) THEN
             atmos%c(i,j)=val
             atmos%c(j,i)= - val
          ENDIF
       END DO
    END DO
  END SUBROUTINE calculate_c_atm

  !> Forcing of the ocean on the atmosphere.
  !> 
  !> \f$ d_{i,j} = (F_i, \nabla^2 \eta_j)\f$ .
  !
  !> @remark
  !> Atmospheric s tensor and oceanic M tensor must be computed before
  !> calling this routine !
  SUBROUTINE calculate_d
    INTEGER :: i,j
    INTEGER :: AllocStat 

    IF ((.NOT. ALLOCATED(atmos%s)) .OR. (.NOT. ALLOCATED(ocean%M))) THEN
       STOP "*** atmos%s and ocean%M must be defined before calling calculate_d ! ***"
    END IF


    IF (natm == 0 ) THEN
       STOP "*** Problem with calculate_d : natm==0 ! ***"
    ELSE
       IF (.NOT. ALLOCATED(atmos%d)) THEN
          ALLOCATE(atmos%d(natm,noc), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       END IF
    END IF
    atmos%d=0.D0

    DO i=1,natm
       DO j=1,noc
          atmos%d(i,j)=atmos%s(i,j) * ocean%M(j,j)
       END DO
    END DO
  END SUBROUTINE calculate_d

  !> Temperature advection terms (atmospheric)
  !> 
  !> \f$ g_{i,j,k} = (F_i, J(F_j, F_k))\f$ .
  !>
  !> and Streamfunction advection terms (atmospheric)
  !> 
  !> \f$ b_{i,j,k} = (F_i, J(F_j, \nabla^2 F_k))\f$ .
  !>
  !> @remark
  !> This is a strict function: it only accepts AKL KKL and LLL types.
  !> For any other combination, it will not calculate anything.
  SUBROUTINE calculate_bg
    INTEGER :: i,j,k
    TYPE(atm_wavenum) :: Ti, Tj, Tk
    REAL(KIND=8) :: val,vb1, vb2, vs1, vs2, vs3, vs4
    INTEGER :: AllocStat
    INTEGER, DIMENSION(3) :: a,b
    INTEGER, DIMENSION(3,3) :: w
    CHARACTER, DIMENSION(3) :: s
    INTEGER :: par,l

    IF (natm == 0 ) THEN
       STOP "*** Problem with calculate_bg : natm==0 ! ***"
    ELSE
       IF (.NOT. ALLOCATED(ipbuf)) THEN
          STOP "*** Problem with calculate_bg : ipbuf not allocated ! ***"
       END IF
       IF (.NOT. ALLOCATED(atmos%g)) THEN
          ALLOCATE(atmos%g(ndim), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       END IF
    END IF

    DO i=1,natm
       DO j=i,natm
          DO k=j,natm

             Ti = awavenum(i)
             Tj = awavenum(j)
             Tk = awavenum(k)

             a(1)=i
             a(2)=j
             a(3)=k

             val=0.D0

             IF ((Ti%typ == "L") .AND. (Tj%typ == "L") .AND. (Tk%typ == "L")) THEN

                CALL piksrt(3,a,par)

                Ti = awavenum(a(1))
                Tj = awavenum(a(2))
                Tk = awavenum(a(3))

                vs3 = S3(Tj%P,Tk%P,Tj%H,Tk%H)
                vs4 = S4(Tj%P,Tk%P,Tj%H,Tk%H)
                val = vs3 * ((delta(Tk%H - Tj%H - Ti%H) - delta(Tk%H &
                     &- Tj%H + Ti%H)) * delta(Tk%P + Tj%P - Ti%P) +&
                     & delta(Tk%H + Tj%H - Ti%H) * (delta(Tk%P - Tj%P&
                     & + Ti%P) - delta(Tk%P - Tj%P - Ti%P))) + vs4 *&
                     & ((delta(Tk%H + Tj%H - Ti%H) * delta(Tk%P - Tj&
                     &%P - Ti%P)) + (delta(Tk%H - Tj%H + Ti%H) -&
                     & delta(Tk%H - Tj%H - Ti%H)) * (delta(Tk%P - Tj&
                     &%P - Ti%P) - delta(Tk%P - Tj%P + Ti%P)))
             ELSE

                s(1)=Ti%typ
                s(2)=Tj%typ
                s(3)=Tk%typ

                w(1,:)=isin("A",s)
                w(2,:)=isin("K",s)
                w(3,:)=isin("L",s)

                IF (ANY(w(1,:)/=0) .AND. ANY(w(2,:)/=0) .AND. ANY(w(3,:)/=0)) THEN
                   b=w(:,1)
                   Ti = awavenum(a(b(1)))
                   Tj = awavenum(a(b(2)))
                   Tk = awavenum(a(b(3)))
                   call piksrt(3,b,par)
                   vb1 = B1(Ti%P,Tj%P,Tk%P)
                   vb2 = B2(Ti%P,Tj%P,Tk%P)
                   val = -2 * sqrt(2.) / pi * Tj%M * delta(Tj%M - Tk%H) * flambda(Ti%P + Tj%P + Tk%P)
                   IF (val /= 0.D0) val = val * (vb1**2 / (vb1**2 - 1) - vb2**2 / (vb2**2 - 1))
                ELSEIF ((w(2,2)/=0) .AND. (w(2,3)==0) .AND. ANY(w(3,:)/=0)) THEN
                   Ti = awavenum(a(w(2,1)))
                   Tj = awavenum(a(w(2,2)))
                   Tk = awavenum(a(w(3,1)))
                   b(1)=w(2,1)
                   b(2)=w(2,2)
                   b(3)=w(3,1)
                   call piksrt(3,b,par)
                   vs1 = S1(Tj%P,Tk%P,Tj%M,Tk%H)
                   vs2 = S2(Tj%P,Tk%P,Tj%M,Tk%H)
                   val = vs1 * (delta(Ti%M - Tk%H - Tj%M) * delta(Ti%P -&
                        & Tk%P + Tj%P) - delta(Ti%M- Tk%H - Tj%M) *&
                        & delta(Ti%P + Tk%P - Tj%P) + (delta(Tk%H - Tj%M&
                        & + Ti%M) + delta(Tk%H - Tj%M - Ti%M)) *&
                        & delta(Tk%P + Tj%P - Ti%P)) + vs2 * (delta(Ti%M&
                        & - Tk%H - Tj%M) * delta(Ti%P - Tk%P - Tj%P) +&
                        & (delta(Tk%H - Tj%M - Ti%M) + delta(Ti%M + Tk%H&
                        & - Tj%M)) * (delta(Ti%P - Tk%P + Tj%P) -&
                        & delta(Tk%P - Tj%P + Ti%P)))
                ENDIF
             ENDIF
             val=par*val*n
             IF (val /= 0.D0) THEN
                CALL add_check(ipbuf,i,j,k,val,atmos%g)
                CALL add_check(ipbuf,j,k,i,val,atmos%g)
                CALL add_check(ipbuf,k,i,j,val,atmos%g)
                CALL add_check(ipbuf,i,k,j,-val,atmos%g)
                CALL add_check(ipbuf,j,i,k,-val,atmos%g)
                CALL add_check(ipbuf,k,j,i,-val,atmos%g)
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    CALL add_to_tensor(ipbuf,atmos%g)
    DO l=1,natm
       ipbuf(l)%nelems=0
    ENDDO
  END SUBROUTINE calculate_bg

  !> Forcing (thermal) of the ocean on the atmosphere.
  !> 
  !> \f$ s_{i,j} = (F_i, \eta_j)\f$ .
  SUBROUTINE calculate_s
    INTEGER :: i,j
    TYPE(atm_wavenum) :: Ti
    TYPE(ocean_wavenum) :: Dj
    REAL(KIND=8) :: val
    INTEGER :: AllocStat 
    IF (natm == 0 ) THEN
       STOP "*** Problem with calculate_s : natm==0 ! ***"
    ELSEIF (noc == 0) then
       STOP "*** Problem with calculate_s : noc==0 ! ***"
    ELSE
       IF (.NOT. ALLOCATED(atmos%s)) THEN
          ALLOCATE(atmos%s(natm,noc), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       END IF
    END IF
    atmos%s=0.D0

    DO i=1,natm
       DO j=1,noc
          Ti = awavenum(i)
          Dj = owavenum(j)
          val=0.D0
          IF (Ti%typ == "A") THEN
             val = flambda(Dj%H) * flambda(Dj%P + Ti%P)
             IF (val /= 0.D0) THEN
                val = val*8*sqrt(2.)*Dj%P/(pi**2 * (Dj%P**2 - Ti%P**2) * Dj%H)
             END IF
          ELSEIF (Ti%typ == "K") THEN
             val = flambda(2 * Ti%M + Dj%H) * delta(Dj%P - Ti%P)
             IF (val /= 0.D0) THEN
                val = val*4*Dj%H/(pi * (-4 * Ti%M**2 + Dj%H**2))
             END IF
          ELSEIF (Ti%typ == "L") THEN
             val = delta(Dj%P - Ti%P) * delta(2 * Ti%H - Dj%H)
          END IF
          IF (val /= 0.D0) THEN
             atmos%s(i,j)=val
          ENDIF
       END DO
    END DO
  END SUBROUTINE calculate_s

  !--------------------------------------------------------!
  ! 2. Inner products in the equations for the ocean       !
  !--------------------------------------------------------!

  !> Forcing of the atmosphere on the ocean.
  !> 
  !> \f$ K_{i,j} = (\eta_i, \nabla^2 F_j)\f$ .
  !
  !> @remark
  !> atmospheric a and s tensors must be computed before calling
  !> this function !
  SUBROUTINE calculate_K
    INTEGER :: i,j
    INTEGER :: AllocStat 

    IF ((.NOT. ALLOCATED(atmos%a)) .OR. (.NOT. ALLOCATED(atmos%s))) THEN
       STOP "*** atmos%a and atmos%s must be defined before calling calculate_K ! ***"
    END IF

    IF (noc == 0 ) THEN
       STOP "*** Problem with calculate_K : noc==0 ! ***"
    ELSEIF (natm == 0 ) THEN
       STOP "*** Problem with calculate_K : natm==0 ! ***"
    ELSE
       IF (.NOT. ALLOCATED(ocean%K)) THEN
          ALLOCATE(ocean%K(noc,natm), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       END IF
    END IF
    ocean%K=0.D0

    DO i=1,noc
       DO j=1,natm
          ocean%K(i,j) = atmos%s(j,i) * atmos%a(j,j)
       END DO
    END DO
  END SUBROUTINE calculate_K

  !> Forcing of the ocean fields on the ocean.
  !> 
  !> \f$ M_{i,j} = (eta_i, \nabla^2 \eta_j)\f$ .
  SUBROUTINE calculate_M
    INTEGER :: i
    TYPE(ocean_wavenum) :: Di
    INTEGER :: AllocStat 
    IF (noc == 0 ) THEN
       STOP "*** Problem with calculate_M : noc==0 ! ***"
    ELSE
       IF (.NOT. ALLOCATED(ocean%M)) THEN
          ALLOCATE(ocean%M(noc,noc), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       END IF
    END IF
    ocean%M=0.D0

    DO i=1,noc
       Di = owavenum(i)
       ocean%M(i,i) = -(n**2) * Di%Nx**2 - Di%Ny**2
    END DO
  END SUBROUTINE calculate_M

  !> Beta term for the ocean
  !> 
  !> \f$ N_{i,j} = (\eta_i, \partial_x \eta_j) \f$.
  SUBROUTINE calculate_N
    INTEGER :: i,j
    TYPE(ocean_wavenum) :: Di,Dj
    REAL(KIND=8) :: val
    INTEGER :: AllocStat 
    IF (noc == 0 ) THEN
       STOP "*** Problem with calculate_N : noc==0 ! ***"
    ELSE
       IF (.NOT. ALLOCATED(ocean%N)) THEN
          ALLOCATE(ocean%N(noc,noc), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       END IF
    END IF
    ocean%N=0.D0

    DO i=1,noc
       DO j=i,noc
          Di = owavenum(i)
          Dj = owavenum(j)
          val = delta(Di%P - Dj%P) * flambda(Di%H + Dj%H)
          IF ((val /= 0.D0).AND.(Di%H/=Dj%H)) THEN
             ocean%N(i,j) = val * (-2) * Dj%H * Di%H * n / ((Dj%H**2 - Di%H**2) * pi)
             ocean%N(j,i) = -val * (-2) * Dj%H * Di%H * n / ((Dj%H**2 - Di%H**2) * pi)
          ENDIF
       END DO
    END DO
  END SUBROUTINE calculate_N

  !> Temperature advection term (passive scalar)
  !> 
  !> \f$ O_{i,j,k} = (\eta_i, J(\eta_j, \eta_k))\f$ .
  !>
  !> and Streamfunction advection terms (oceanic)
  !> 
  !> \f$ C_{i,j,k} = (\eta_i, J(\eta_j,\nabla^2 \eta_k))\f$ .
  SUBROUTINE calculate_OC
    INTEGER :: i,j,k
    REAL(KIND=8) :: vs3,vs4,val
    TYPE(ocean_wavenum) :: Di,Dj,Dk
    INTEGER :: AllocStat
    INTEGER :: l
    IF (noc == 0 ) THEN
       STOP "*** Problem with calculate_O : noc==0 ! ***"
    ELSE
       IF (.NOT. ALLOCATED(ipbuf)) THEN
          STOP "*** Problem with calculate_OC : ipbuf not allocated ! ***"
       END IF
       IF (.NOT. ALLOCATED(ocean%O)) THEN
          ALLOCATE(ocean%O(ndim), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       END IF
    END IF

    DO i=1,noc
       DO j=i,noc
          DO k=j,noc
             Di = owavenum(i)
             Dj = owavenum(j)
             Dk = owavenum(k)
             vs3 = S3(Dj%P,Dk%P,Dj%H,Dk%H)
             vs4 = S4(Dj%P,Dk%P,Dj%H,Dk%H)
             val = vs3*((delta(Dk%H - Dj%H - Di%H) - delta(Dk%H - Dj&
                  &%H + Di%H)) * delta(Dk%P + Dj%P - Di%P) + delta(Dk&
                  &%H + Dj%H - Di%H) * (delta(Dk%P - Dj%P + Di%P) -&
                  & delta(Dk%P - Dj%P - Di%P))) + vs4 * ((delta(Dk%H &
                  &+ Dj%H - Di%H) * delta(Dk%P - Dj%P - Di%P)) +&
                  & (delta(Dk%H - Dj%H + Di%H) - delta(Dk%H - Dj%H -&
                  & Di%H)) * (delta(Dk%P - Dj%P - Di%P) - delta(Dk%P &
                  &- Dj%P + Di%P)))
             val = val * n / 2
             IF (val /= 0.D0) THEN
                CALL add_check(ipbuf,i,j,k,val,ocean%O) 
                CALL add_check(ipbuf,j,k,i,val,ocean%O)
                CALL add_check(ipbuf,k,i,j,val,ocean%O)
                CALL add_check(ipbuf,i,k,j,-val,ocean%O)
                CALL add_check(ipbuf,j,i,k,-val,ocean%O)
                CALL add_check(ipbuf,k,j,i,-val,ocean%O)
             END IF
          END DO
       END DO
    END DO

    CALL add_to_tensor(ipbuf,ocean%O)    
    DO l=1,noc
       ipbuf(l)%nelems=0
    ENDDO
  END SUBROUTINE calculate_OC


  !> Short-wave radiative forcing of the ocean
  !> 
  !> \f$ W_{i,j} = (\eta_i, F_j)\f$ .
  !
  !> @remark
  !> atmospheric s tensor must be computed before calling
  !> this function !
  SUBROUTINE calculate_W
    INTEGER :: i,j
    INTEGER :: AllocStat 

    IF (.NOT. ALLOCATED(atmos%s)) THEN
       STOP "*** atmos%s must be defined before calling calculate_W ! ***"
    END IF

    IF (noc == 0 ) THEN
       STOP "*** Problem with calculate_W : noc==0 ! ***"
    ELSEIF (natm == 0 ) THEN
       STOP "*** Problem with calculate_W : natm==0 ! ***"
    ELSE
       IF (.NOT. ALLOCATED(ocean%W)) THEN
          ALLOCATE(ocean%W(noc,natm), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       END IF
    END IF
    ocean%W=0.D0

    DO i=1,noc
       DO j=1,natm
          ocean%W(i,j) = atmos%s(j,i)
       END DO
    END DO
  END SUBROUTINE calculate_W

  !-----------------------------------------------------!
  !                                                     !
  ! Initialisation routine                              !
  !                                                     !
  !-----------------------------------------------------!

  !> Initialisation of the inner product
  SUBROUTINE init_inprod
    INTEGER :: i,j
    INTEGER :: AllocStat
    LOGICAL :: ex

    ! Definition of the types and wave numbers tables

    ALLOCATE(owavenum(noc),awavenum(natm), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    j=0
    DO i=1,nbatm
       IF (ams(i,1)==1) THEN
          awavenum(j+1)%typ='A'
          awavenum(j+2)%typ='K'
          awavenum(j+3)%typ='L'

          awavenum(j+1)%P=ams(i,2)
          awavenum(j+2)%M=ams(i,1)
          awavenum(j+2)%P=ams(i,2)
          awavenum(j+3)%H=ams(i,1)
          awavenum(j+3)%P=ams(i,2)

          awavenum(j+1)%Ny=REAL(ams(i,2))
          awavenum(j+2)%Nx=REAL(ams(i,1))
          awavenum(j+2)%Ny=REAL(ams(i,2))
          awavenum(j+3)%Nx=REAL(ams(i,1))
          awavenum(j+3)%Ny=REAL(ams(i,2))

          j=j+3
       ELSE
          awavenum(j+1)%typ='K'
          awavenum(j+2)%typ='L'

          awavenum(j+1)%M=ams(i,1)
          awavenum(j+1)%P=ams(i,2)
          awavenum(j+2)%H=ams(i,1)
          awavenum(j+2)%P=ams(i,2)

          awavenum(j+1)%Nx=REAL(ams(i,1))
          awavenum(j+1)%Ny=REAL(ams(i,2))
          awavenum(j+2)%Nx=REAL(ams(i,1))
          awavenum(j+2)%Ny=REAL(ams(i,2))

          j=j+2

       ENDIF
    ENDDO

    DO i=1,noc
       owavenum(i)%H=oms(i,1)
       owavenum(i)%P=oms(i,2)

       owavenum(i)%Nx=oms(i,1)/2.D0
       owavenum(i)%Ny=oms(i,2)

    ENDDO

    ! Computation of the atmospheric inner products tensors

    ! Allocating the buffer
    ALLOCATE(ipbuf(ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    DO i=1,ndim
       ALLOCATE(ipbuf(i)%elems(1000), STAT=AllocStat)
       IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    END DO

    CALL calculate_a
    INQUIRE(FILE='atmos_g.ipf',EXIST=ex)
    IF (ex) THEN
       IF (.NOT. ALLOCATED(atmos%g)) THEN
          ALLOCATE(atmos%g(ndim), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       END IF
       
       CALL load_tensor_from_file('atmos_g.ipf',atmos%g)
    ELSE
       CALL calculate_bg
       CALL write_tensor_to_file('atmos_g.ipf',atmos%g)
    ENDIF

    CALL calculate_s
    CALL calculate_c_atm

    ! Computation of the oceanic inner products tensors

    CALL calculate_M
    CALL calculate_N
    
    INQUIRE(FILE='ocean_O.ipf',EXIST=ex)
    IF (ex) THEN
       IF (.NOT. ALLOCATED(ocean%O)) THEN
          ALLOCATE(ocean%O(ndim), STAT=AllocStat)
          IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
       END IF

       CALL load_tensor_from_file('ocean_O.ipf',ocean%O)
    ELSE
       CALL calculate_OC
       CALL write_tensor_to_file('ocean_O.ipf',ocean%O)
    ENDIF


    CALL calculate_W
    CALL calculate_K

    ! A last atmospheric one that needs ocean%M

    CALL calculate_d

    DEALLOCATE(ipbuf, STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Deallocation problem ! ***"

  END SUBROUTINE init_inprod


END MODULE inprod_analytic

