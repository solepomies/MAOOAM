
!  maooam_lyap.f90
!
!> Fortran 90 implementation of the modular arbitrary-order ocean-atmosphere
!> model MAOOAM computing the Lyapunov spectrum.
!
!> @copyright                                                               
!> 2016 Lesley De Cruz, Sebastian Schubert & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM maooam_lyap
  USE params, only: ndim, dt, tw, tw_snap, t_trans, t_run, writeout, rescaling_time
  USE aotensor_def, only: init_aotensor
  USE tl_ad_tensor, only: init_tltensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE tl_ad_integrator, only: init_tl_ad_integrator,prop_step
  USE lyap_vectors, only: lyapunov,loclyap,init_lyap,multiply_prop,benettin_step
  USE stat
  USE lyap_stat
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X          !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xwrite     !< Copy of the latter to be written
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew       !< Updated state variable
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: prop_buf !< Buffer for the propagator
  REAL(KIND=8) :: t=0.D0                                !< Time variable
  REAL(KIND=8) :: t_up
  INTEGER :: IndexBen,WRSTAT
  CHARACTER(LEN=19) :: FMTX
  INTEGER :: i, next, IndexSnap, WRSTAT2
  CHARACTER(LEN=9) :: arg
  LOGICAL :: cont_evol    !< True if the initial state is to be read in snapshot_trans.dat (i.e. the previous evolution is to be continued)
  LOGICAL :: ex

  PRINT*, 'Model MAOOAM v1.0'
  PRINT*, '      - with computation of the Lyapunov spectrum'
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

  CALL init_aotensor    ! Compute the tensors
  CALL init_tltensor   
  CALL load_IC          ! Load the initial condition

  CALL init_integrator        ! Initialize the integrator
  CALL init_tl_ad_integrator  ! Initialize tangent linear integrator
  CALL init_lyap              ! Initialize Lyapunov computation
  write(FMTX,'(A10,i3,A6)') '(F10.2,4x,',ndim,'E15.5)'
  t_up=dt/t_trans*100.D0

  IF (writeout) THEN
     OPEN(10,file='evol_field.dat')
     OPEN(11,file='lyapunov_exponents.dat',status='replace',form='UNFORMATTED',access='DIRECT',recl=8*ndim)
  END IF

  ALLOCATE(X(0:ndim),Xnew(0:ndim),prop_buf(ndim,ndim))
  X=IC
  IF (cont_evol) THEN
     INQUIRE(FILE='snapshots_trans.dat',EXIST=ex,NEXTREC=next)
     IF (ex) THEN
        OPEN(12,file='snapshots_trans.dat')
        READ(12,REC=next-1) X
        CLOSE(12)
     END IF
  END IF

  IF (writeout) OPEN(12,file='snapshots_trans.dat',status='replace',form='UNFORMATTED',access='DIRECT',recl=8*ndim)


  PRINT*, 'Starting the transient time evolution... t_trans = ',t_trans

  IndexSnap=0
  DO WHILE (t<t_trans)
     CALL step(X,t,dt,Xnew)
     X=Xnew
     IF (mod(t/t_trans*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_trans*100.D0,char(13)
     IF (writeout .AND. mod(t,tw_snap)<dt) THEN
        IndexSnap=IndexSnap+1
        WRITE(12,rec=IndexSnap,iostat=WRSTAT2) X
     END IF
  END DO

  IF (writeout) CLOSE(12)

  PRINT*, 'Starting the time evolution... t_run = ',t_run

  CALL init_stat
  CALL lyap_init_stat
  t=0.D0
  IndexBen=0
  t_up=dt/t_run*100.D0

  DO WHILE (t<t_run)
     CALL prop_step(X,prop_buf,t,dt,Xnew,.false.) ! Obtains propagator prop_buf at X
     CALL multiply_prop(prop_buf) ! Multiplies prop_buf with prop
     X=Xnew
     IF (mod(t,rescaling_time)<dt) THEN
        CALL  benettin_step ! Performs QR step with prop
        CALL lyap_acc(loclyap)
     END IF
     IF (mod(t,tw)<dt) THEN
        CALL acc(X)
        Xwrite=X
        DO i=1,ndim
           IF (abs(Xwrite(i))<1.D-50) Xwrite(i)=0
        END DO
        IF (writeout) WRITE(10,FMTX) t,Xwrite(1:ndim)
        IndexBen=IndexBen+1
        IF (writeout) WRITE(11,rec=IndexBen,iostat=WRSTAT) loclyap
     END IF
     IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
  END DO
  PRINT*, 'Evolution finished.'

  IF (writeout) CLOSE(10)
  IF (writeout) CLOSE(11)

  IF (writeout) THEN
     OPEN(10,file='mean_lyapunov.dat')
     lyapunov=lyap_mean()
     WRITE(10,*) 'mean',lyapunov(1:ndim)
     lyapunov=lyap_var()
     WRITE(10,*) 'var',lyapunov(1:ndim)
  END IF

  IF (writeout) THEN
     OPEN(10,file='mean_field.dat')
     X=mean()
     WRITE(10,*) 'mean',X(1:ndim)
     X=var()
     WRITE(10,*) 'var',X(1:ndim)
  END IF

END PROGRAM maooam_lyap
