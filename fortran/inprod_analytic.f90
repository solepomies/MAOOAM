
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

  USE params, only: nbatm, nboc, nball, natm, noc, nall, n, oms, ams, allms, pi, Hamax, Pamax, Homax, Pomax
  USE util, only: isin,piksrt
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_inprod

  !> General block specification type
  TYPE :: all_wavenum
     CHARACTER :: typ
     INTEGER :: M=0,P=0,H=0
     REAL(KIND=8) :: Nx=0.,Ny=0.
  END TYPE all_wavenum

  !> Type holding the atmospheric inner products tensors
  TYPE :: atm_tensors
     PROCEDURE(calculate_a), POINTER, NOPASS :: a
     PROCEDURE(calculate_b), POINTER, NOPASS :: b
     PROCEDURE(calculate_c_atm), POINTER, NOPASS :: c
     PROCEDURE(calculate_g), POINTER, NOPASS :: g
  END TYPE atm_tensors

  !> Atmospheric blocks specification
  TYPE(all_wavenum), DIMENSION(:), ALLOCATABLE, PUBLIC :: awavenum 
  !> Oceanic blocks specification
  TYPE(all_wavenum), DIMENSION(:), ALLOCATABLE, PUBLIC :: owavenum 
  !> All blocks specification
  TYPE(all_wavenum), DIMENSION(:), ALLOCATABLE, PUBLIC :: allwavenum 

  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: correspatm  !< This list will give where to find the i-th atmospheric basis function in the list of all basis functions
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: correspoc   !< This list will give where to find the i-th oceanic basis function in the list of all basis functions

  !> Atmospheric tensors
  TYPE(atm_tensors), PUBLIC :: atmos


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
  !-----------------------------------------------------!
  
  !> Eigenvalues of the Laplacian (atmospheric)
  !> 
  !> \f$ a_{i,j} = (F_i, \nabla^2 F_j)\f$ .
  REAL(KIND=8) FUNCTION calculate_a(i,j)
    INTEGER, INTENT(IN) :: i,j
    TYPE(all_wavenum) :: Ti
    
    calculate_a = 0.D0
    IF (i==j) THEN
       Ti = allwavenum(i)
       calculate_a = -(n**2) * Ti%Nx**2 - Ti%Ny**2
    ENDIF
  END FUNCTION calculate_a

  !> Streamfunction advection terms (atmospheric)
  !> 
  !> \f$ b_{i,j,k} = (F_i, J(F_j, \nabla^2 F_k))\f$ .
  REAL(KIND=8) FUNCTION calculate_b(i,j,k)
    INTEGER, INTENT(IN) :: i,j,k

    calculate_b = calculate_a(k,k) * calculate_g(i,j,k)

  END FUNCTION calculate_b

  !> Beta term for the atmosphere
  !> 
  !> \f$ c_{i,j} = (F_i, \partial_x F_j)\f$ .
  REAL(KIND=8) FUNCTION calculate_c_atm(i,j)
    INTEGER, INTENT(IN) :: i,j
    TYPE(all_wavenum) :: Ti, Tj

    Ti = allwavenum(i)
    Tj = allwavenum(j)
    calculate_c_atm = 0.D0
    IF ((Ti%typ == "K") .AND. (Tj%typ == "L")) THEN 
       calculate_c_atm = n * Ti%M * delta(Ti%M - Tj%H) * delta(Ti%P - Tj%P)
    ELSE IF ((Ti%typ == "L") .AND. (Tj%typ == "K")) THEN
       Ti = allwavenum(j)
       Tj = allwavenum(i)
       calculate_c_atm = - n * Ti%M * delta(Ti%M - Tj%H) * delta(Ti%P - Tj%P)
    END IF
    
  END FUNCTION calculate_c_atm

  !> Temperature advection terms (atmospheric)
  !> 
  !> \f$ g_{i,j,k} = (F_i, J(F_j, F_k))\f$ .
  REAL(KIND=8) FUNCTION calculate_g(i,j,k)
    INTEGER, INTENT(IN) :: i,j,k
    TYPE(all_wavenum) :: Ti,Tj,Tk
    REAL(KIND=8) :: val,vb1, vb2, vs1, vs2, vs3, vs4
    INTEGER, DIMENSION(3) :: a,b
    INTEGER, DIMENSION(3,3) :: w
    CHARACTER, DIMENSION(3) :: s
    INTEGER :: par

    Ti = allwavenum(i)
    Tj = allwavenum(j)
    Tk = allwavenum(k)

    a(1)=i
    a(2)=j
    a(3)=k

    par=0.D0
    val=0.D0

    IF ((Ti%typ == "L") .AND. (Tj%typ == "L") .AND. (Tk%typ == "L")) THEN
       
       CALL piksrt(3,a,par)

       Ti = allwavenum(a(1))
       Tj = allwavenum(a(2))
       Tk = allwavenum(a(3))

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
          Ti = allwavenum(a(b(1)))
          Tj = allwavenum(a(b(2)))
          Tk = allwavenum(a(b(3)))
          call piksrt(3,b,par)
          vb1 = B1(Ti%P,Tj%P,Tk%P)
          vb2 = B2(Ti%P,Tj%P,Tk%P)
          val = -2 * sqrt(2.) / pi * Tj%M * delta(Tj%M - Tk%H) * flambda(Ti%P + Tj%P + Tk%P)
          IF (val /= 0.D0) val = val * (vb1**2 / (vb1**2 - 1) - vb2**2 / (vb2**2 - 1))
       ELSEIF ((w(2,2)/=0) .AND. (w(2,3)==0) .AND. ANY(w(3,:)/=0)) THEN
          Ti = allwavenum(a(w(2,1)))
          Tj = allwavenum(a(w(2,2)))
          Tk = allwavenum(a(w(3,1)))
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
    calculate_g=par*val*n
 
  END FUNCTION calculate_g

  !-----------------------------------------------------!
  !                                                     !
  ! Initialisation routine                              !
  !                                                     !
  !-----------------------------------------------------!

  !> Initialisation of the inner product
  SUBROUTINE init_inprod
    INTEGER :: i,j
    INTEGER :: AllocStat
    INTEGER :: H,P

    IF (natm == 0 ) THEN
       STOP "*** Problem : natm==0 ! ***"
    ELSEIF (noc == 0) then
       STOP "*** Problem : noc==0 ! ***"
    END IF


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

    j=0
    DO i=1,nboc
       IF (oms(i,1)==1) THEN
          owavenum(j+1)%typ='A'
          owavenum(j+2)%typ='K'
          owavenum(j+3)%typ='L'

          owavenum(j+1)%P=oms(i,2)
          owavenum(j+2)%M=oms(i,1)
          owavenum(j+2)%P=oms(i,2)
          owavenum(j+3)%H=oms(i,1)
          owavenum(j+3)%P=oms(i,2)

          owavenum(j+1)%Ny=REAL(oms(i,2))
          owavenum(j+2)%Nx=REAL(oms(i,1))
          owavenum(j+2)%Ny=REAL(oms(i,2))
          owavenum(j+3)%Nx=REAL(oms(i,1))
          owavenum(j+3)%Ny=REAL(oms(i,2))

          j=j+3
       ELSE
          owavenum(j+1)%typ='K'
          owavenum(j+2)%typ='L'

          owavenum(j+1)%M=oms(i,1)
          owavenum(j+1)%P=oms(i,2)
          owavenum(j+2)%H=oms(i,1)
          owavenum(j+2)%P=oms(i,2)

          owavenum(j+1)%Nx=REAL(oms(i,1))
          owavenum(j+1)%Ny=REAL(oms(i,2))
          owavenum(j+2)%Nx=REAL(oms(i,1))
          owavenum(j+2)%Ny=REAL(oms(i,2))

          j=j+2

       ENDIF

    ENDDO

    ALLOCATE(allwavenum(nall), correspatm(natm), correspoc(noc), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    j=0
    DO i=1,nball
       H=allms(i,1)
       P=allms(i,2)
       IF (H==1) THEN
          allwavenum(j+1)%typ='A'
          allwavenum(j+2)%typ='K'
          allwavenum(j+3)%typ='L'

          allwavenum(j+1)%P=P
          allwavenum(j+2)%M=H
          allwavenum(j+2)%P=P
          allwavenum(j+3)%H=H
          allwavenum(j+3)%P=P

          allwavenum(j+1)%Ny=REAL(P)
          allwavenum(j+2)%Nx=REAL(H)
          allwavenum(j+2)%Ny=REAL(P)
          allwavenum(j+3)%Nx=REAL(H)
          allwavenum(j+3)%Ny=REAL(P)

          IF ((H .le. Hamax) .AND. (P .le. Pamax)) THEN
             correspatm((P-1)*3+1)=j+1
             correspatm((P-1)*3+2)=j+2
             correspatm((P-1)*3+3)=j+3
          ENDIF

          IF ((H .le. Homax) .AND. (P .le. Pomax)) THEN
             correspoc((P-1)*3+1)=j+1
             correspoc((P-1)*3+2)=j+2
             correspoc((P-1)*3+3)=j+3
          ENDIF

          j=j+3
       ELSE
          allwavenum(j+1)%typ='K'
          allwavenum(j+2)%typ='L'

          allwavenum(j+1)%M=H
          allwavenum(j+1)%P=P
          allwavenum(j+2)%H=H
          allwavenum(j+2)%P=P

          allwavenum(j+1)%Nx=REAL(H)
          allwavenum(j+1)%Ny=REAL(P)
          allwavenum(j+2)%Nx=REAL(H)
          allwavenum(j+2)%Ny=REAL(P)

          IF ((H .le. Hamax) .AND. (P .le. Pamax)) THEN
             correspatm((3+2*(H-2))*Pamax+(P-1)*2+1)=j+1
             correspatm((3+2*(H-2))*Pamax+(P-1)*2+2)=j+2
          ENDIF

          IF ((H .le. Homax) .AND. (P .le. Pomax)) THEN
             correspoc((3+2*(H-2))*Pomax+(P-1)*2+1)=j+1
             correspoc((3+2*(H-2))*Pomax+(P-1)*2+2)=j+2
          ENDIF

          j=j+2

       ENDIF

    ENDDO

    ! Pointing to the inner products functions

    atmos%a => calculate_a
    atmos%g => calculate_g
    atmos%b => calculate_b
    atmos%c => calculate_c_atm

  END SUBROUTINE init_inprod


END MODULE inprod_analytic

