
!  maooam.f90
!
!> Fortran 90 implementation of the modular arbitrary-order ocean-atmosphere 
!> model MAOOAM.
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM maooam 
  USE params, only: ndim, dt, tw, tw_snap, t_trans, t_run, writeout
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE stat
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X       !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew    !< Updated state variable
  REAL(KIND=8) :: t=0.D0                             !< Time variable
  REAL(KIND=8) :: t_up
  INTEGER :: i, next, IndexSnap, WRSTAT
  CHARACTER(LEN=9) :: arg
  LOGICAL :: cont_evol    !< True if the initial state is to be read in snapshot_trans.dat (i.e. the previous evolution is to be continued)
  LOGICAL :: ex

  PRINT*, 'Model MAOOAM v1.3'
  PRINT*, 'Loading information...'

  ! Initializing cont_evol
  cont_evol=.FALSE.
  arg='arg'
  i=0
  DO WHILE ((.NOT. cont_evol) .AND. LEN_TRIM(arg) /= 0)
     CALL get_command_argument(i, arg)
     IF (TRIM(arg)=='continue') cont_evol=.TRUE.
     i=i+1
  END DO

  CALL init_aotensor    ! Compute the tensor

  CALL load_IC          ! Load the initial condition

  CALL init_integrator  ! Initialize the integrator

  t_up=dt/t_trans*100.D0

  IF (writeout) OPEN(10,file='evol_field.dat')

  ALLOCATE(X(0:ndim),Xnew(0:ndim))

  X=IC
  IF (cont_evol) THEN
     INQUIRE(FILE='snapshots_trans.dat',EXIST=ex,NEXTREC=next)
     IF (ex) THEN
        OPEN(11,file='snapshots_trans.dat')
        READ(11,REC=next-1) X
        CLOSE(11)
     END IF
  END IF

  IF (writeout) OPEN(11,file='snapshots_trans.dat',status='replace',form='UNFORMATTED',access='DIRECT',recl=8*ndim)

  PRINT*, 'Starting the transient time evolution...'

  IndexSnap=0
  DO WHILE (t<t_trans)
     CALL step(X,t,dt,Xnew)
     X=Xnew
     IF (mod(t/t_trans*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_trans*100.D0,char(13)
     IF (writeout .AND. mod(t,tw_snap)<dt) THEN
        IndexSnap=IndexSnap+1
        WRITE(11,rec=IndexSnap,iostat=WRSTAT) X
     END IF
  END DO

  IF (writeout) CLOSE(11)

  PRINT*, 'Starting the time evolution...'

  CALL init_stat
  
  t=0.D0
  t_up=dt/t_run*100.D0

  IF (writeout) WRITE(10,*) t,X(1:ndim)

  DO WHILE (t<t_run)
     CALL step(X,t,dt,Xnew)
     X=Xnew
     IF (mod(t,tw)<dt) THEN
        IF (writeout) WRITE(10,*) t,X(1:ndim)
        CALL acc(X)
     END IF
     IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
  END DO

  PRINT*, 'Evolution finished.'

  IF (writeout) CLOSE(10)

  IF (writeout) OPEN(10,file='mean_field.dat')

  X=mean()
  IF (writeout) WRITE(10,*) X(1:ndim)
  IF (writeout) CLOSE(10)

END PROGRAM maooam 
