
! icdelta_def.f90
!
!>  Module to load the perturbation initial condition.
!
!> @copyright                                                               
!> 2016 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

MODULE icdelta_def

  USE params, only: natm,noc,ndim
  USE util, only: str,rstr
  USE inprod_analytic, only:awavenum,owavenum
  IMPLICIT NONE

  PRIVATE

  LOGICAL :: exists !< Boolean to test for file existence.
  
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC :: ICdelta !< Initial condition vector

  PUBLIC ::load_ICdelta

CONTAINS

  !> Subroutine to load the initial condition if ICdelta.nml exists.
  !> If it does not, then write ICdelta.nml with random initial condition.
  SUBROUTINE load_ICdelta
    INTEGER :: i,AllocStat
    CHARACTER(len=20) :: fm
    REAL(KIND=8) :: size_of_random_noise
    CHARACTER(LEN=4) :: init_type 
    NAMELIST /IClist/ ICdelta
    NAMELIST /RAND/ init_type,size_of_random_noise 



    fm(1:6)='(F3.1)'
   
    IF (ndim == 0) STOP "*** Number of dimensions is 0! ***"
    ALLOCATE(ICdelta(0:ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    INQUIRE(FILE='./ICdelta.nml',EXIST=exists)

    IF (exists) THEN
       OPEN(8, file="ICdelta.nml", status='OLD', recl=80, delim='APOSTROPHE')
       READ(8,nml=IClist)
       READ(8,nml=RAND)
       SELECT CASE (init_type)
         CASE ('rand')
           CALL random_number(ICdelta)
           ICdelta=2*(ICdelta-0.5)
           ICdelta=ICdelta*size_of_random_noise*10.D0
           ICdelta(0)=1.0d0
           WRITE(6,*) "*** ICdelta.nml namelist written. Starting with random initial condition !***"
         CASE ('zero')
           ICdelta=0
           ICdelta(0)=1.0d0
           WRITE(6,*) "*** ICdelta.nml namelist written. Starting with initial condition in ICdelta.nml !***"
         CASE ('read')
           !nothing has to be done ICdelta has already the right values
           WRITE(6,*) "*** ICdelta.nml namelist written. Starting with initial condition in ICdelta.nml !***"
       END SELECT
       CLOSE(8)
    ELSE
       CALL random_number(ICdelta)
       ICdelta=2*(ICdelta-0.5)
       size_of_random_noise=1.D-3
       ICdelta=ICdelta*size_of_random_noise*10.D0
       ICdelta(0)=1.0d0
       init_type="rand"
       WRITE(6,*) "*** ICdelta.nml namelist written. Starting with 0 as initial condition !***"

    END IF
    OPEN(8, file="ICdelta.nml", status='REPLACE')
    WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
    WRITE(8,'(a)') "! Namelist file :                                                              !"
    WRITE(8,'(a)') "! Initial condition.                                                           !"
    WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
    WRITE(8,*) ""
    WRITE(8,'(a)') "&ICLIST"
    WRITE(8,*) " ! psi variables"
    DO i=1,natm
       WRITE(8,*) " ICdelta("//TRIM(str(i))//") = ",ICdelta(i+natm),"   ! typ= "&
            &//awavenum(i)%typ//", Nx= "//TRIM(rstr(awavenum(i)&
            &%Nx,fm))//", Ny= "//TRIM(rstr(awavenum(i)%Ny,fm))
    END DO
    WRITE(8,*) " ! theta variables"
    DO i=1,natm
       WRITE(8,*) " ICdelta("//TRIM(str(i+natm))//") = ",ICdelta(i+natm),"   ! typ= "&
            &//awavenum(i)%typ//", Nx= "//TRIM(rstr(awavenum(i)&
            &%Nx,fm))//", Ny= "//TRIM(rstr(awavenum(i)%Ny,fm))
    END DO

    WRITE(8,*) " ! A variables"
    DO i=1,noc
       WRITE(8,*) " ICdelta("//TRIM(str(i+2*natm))//") = ",ICdelta(i+2*natm),"   ! Nx&
            &= "//TRIM(rstr(owavenum(i)%Nx,fm))//", Ny= "&
            &//TRIM(rstr(owavenum(i)%Ny,fm))
    END DO
    WRITE(8,*) " ! T variables"
    DO i=1,noc
       WRITE(8,*) " ICdelta("//TRIM(str(i+noc+2*natm))//") = ",ICdelta(i+2*natm+noc),"   &
            &! Nx= "//TRIM(rstr(owavenum(i)%Nx,fm))//", Ny= "&
            &//TRIM(rstr(owavenum(i)%Ny,fm))
    END DO

    WRITE(8,'(a)') "&END"
    WRITE(8,*) ""
    WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
    WRITE(8,'(a)') "! Initialisation type.                                                         !"
    WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
    WRITE(8,'(a)') "! type = 'read': use ICdelta; 'rand': random state; 'zero': zero condition "
    WRITE(8,'(a)') "! The seed is specified in IC.nml"
    WRITE(8,*) ""
    WRITE(8,'(a)') "&RAND"
    WRITE(8,'(a)') "  init_type= '"//init_type//"'" 
    WRITE(8,'(a,d15.7)') "  size_of_random_noise = ",size_of_random_noise
    WRITE(8,'(a)') "&END"
    WRITE(8,*) ""
    CLOSE(8)
  END SUBROUTINE load_ICdelta

END MODULE icdelta_def



