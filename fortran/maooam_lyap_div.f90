
!  maooam_lyap_div.f90
!
!> Fortran 90 implementation of the modular arbitrary-order ocean-atmosphere
!> model MAOOAM computing the first Lyapunov exponent with the divergence
!> method.
!
!> @copyright                                                               
!> 2016 Lesley De Cruz, Sebastian Schubert & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM maooam_lyap_div
  USE params, only: ndim, dt, tw, t_trans, t_run, writeout, rescaling_time
  USE aotensor_def, only: init_aotensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE stat
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X,Xp          !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew,Xnewp    !< Updated state variable
  REAL(KIND=8) :: t=0.D0                                   !< Time variable
  REAL(KIND=8) :: resc=1.D-9                               !< Variable rescaling factor for the divergence method
  REAL(KIND=8) :: t_up
  CHARACTER(LEN=19) :: FMTX

  PRINT*, 'Model MAOOAM v1.0'
  PRINT*, '      - with computation of the first Lyapunov exponent with the divergence method'
  PRINT*, 'Loading information...'

  CALL init_aotensor    ! Compute the tensors
  CALL load_IC          ! Load the initial condition

  CALL init_integrator  ! Initialize the integrator
  write(FMTX,'(A10,i3,A6)') '(F10.2,4x,',ndim,'E15.5)'
  t_up=dt/t_trans*100.D0

  IF (writeout) THEN
     OPEN(10,file='evol_field.dat')
     OPEN(12,file='lyapunov_exponents_div.dat')
  END IF

  ALLOCATE(X(0:ndim),Xnew(0:ndim))
  ALLOCATE(Xp(0:ndim),Xnewp(0:ndim))
  X=IC
  PRINT*, 'Starting the transient time evolution... t_trans = ',t_trans

  DO WHILE (t<t_trans)
     CALL step(X,t,dt,Xnew)
     X=Xnew
     IF (mod(t/t_trans*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_trans*100.D0,char(13)
  END DO

  PRINT*, 'Starting the time evolution... t_run = ',t_run

  CALL init_stat
  Xp(0)=1.0d0
  Xp(1:ndim)=X(1:ndim)+resc/sqrt(dble(ndim))
  t=0.D0
  t_up=dt/t_run*100.D0
  DO WHILE (t<t_run)
     CALL step(X,t,dt,Xnew)
     X=Xnew
     CALL step(Xp,t,dt,Xnewp)
     Xp=Xnewp
     t=t-dt ! Time was incremented one step too much
     IF (mod(t,rescaling_time)<dt) THEN
        write(12,*) log(sqrt(sum((X(1:ndim)-Xp(1:ndim))**2.0d0))/resc)/rescaling_time
        Xp(1:ndim)=X(1:ndim)+(Xp(1:ndim)-X(1:ndim))/sqrt(sum((X(1:ndim)-Xp(1:ndim))**2.0d0))*resc
     END IF
     IF (mod(t,tw)<dt) THEN
        CALL acc(X)
        IF (writeout) WRITE(10,FMTX) t,X(1:ndim) 
        CONTINUE
     END IF
   
     IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
  END DO
  PRINT*, 'Evolution finished.'

  IF (writeout) CLOSE(10)
  IF (writeout) CLOSE(12)

  IF (writeout) THEN
     OPEN(10,file='mean_field.dat')
     X=mean()
     WRITE(10,*) 'mean',X(1:ndim)
     X=var()
     WRITE(10,*) 'var',X(1:ndim)
  END IF

END PROGRAM maooam_lyap_div
