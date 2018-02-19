
! tl_ad_integrator.f90
!
!> Tangent Linear (TL) and Adjoint (AD) model versions of MAOOAM.
!> Integrators module.
!
!> @copyright                                                               
!> 2016 Lesley De Cruz, Jonathan Demaeyer & Sebastian Schubert.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!>  @remark                                                                 
!>  This module actually contains the Heun algorithm routines.
!>  The user can modify it according to its preferred integration scheme.
!>  For higher-order schemes, additional buffers will probably have to be defined.
!                                                                           
!---------------------------------------------------------------------------



MODULE tl_ad_integrator

  USE util, only: init_one
  USE params, only: ndim
  USE tensor, only: sparse_mul3
  USE aotensor_def, only: aotensor

  USE tl_ad_tensor, only: ad,tl,jacobian_mat
  IMPLICIT NONE

  PRIVATE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y1 !< Buffer to hold the intermediate position (Heun algorithm) of the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_f0 !< Buffer to hold tendencies at the initial position of the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_f1 !< Buffer to hold tendencies at the intermediate position of the tangent linear model

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y11 !< Buffer to hold the intermediate position (Heun algorithm) of the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_f00 !< Buffer to hold tendencies at the initial position of the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_f11 !< Buffer to hold tendencies at the intermediate position of the tangent linear model

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j1 !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j2 !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j1h !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j2h !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: one        !< unit matrix 

    
  PUBLIC :: init_tl_ad_integrator, ad_step, tl_step, evolve_ad_step, evolve_tl_step, prop_step

CONTAINS

  !> Routine computing the tendencies of the nonlinear model
  !> @param t Time at which the tendencies have to be computed. Actually not needed for autonomous systems.
  !> @param y Point at which the tendencies have to be computed.
  !> @param res vector to store the result.
  !> @remark Note that it is NOT safe to pass `y` as a result bufer, 
  !> as this operation does multiple passes.
  SUBROUTINE tendencies(t,y,res)
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    CALL sparse_mul3(aotensor, y, y, res)
  END SUBROUTINE tendencies

  !> Routine to initialise the integration buffers.
  SUBROUTINE init_tl_ad_integrator
    INTEGER :: AllocStat
    ALLOCATE(buf_y1(0:ndim),buf_f0(0:ndim),buf_f1(0:ndim),&
         &buf_y11(0:ndim),buf_f00(0:ndim),buf_f11(0:ndim),&
         &buf_j1(ndim,ndim),buf_j2(ndim,ndim),one(ndim,ndim),&
         &buf_j1h(ndim,ndim),buf_j2h(ndim,ndim),STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    CALL init_one(one)
  END SUBROUTINE init_tl_ad_integrator



  !-----------------------------------------------------!
  !                                                     !
  ! Adjoint model integrator                            !
  !                                                     !
  !-----------------------------------------------------!
  
  !> Routine to perform an integration step (Heun algorithm) of the adjoint model. The incremented time is returned.
  !> @param y Initial point.
  !> @param ystar Adjoint model at the point ystar.
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param res Final point after the step.
  SUBROUTINE ad_step(y,ystar,t,dt,res)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y,ystar
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    
    CALL ad(t,ystar,y,buf_f0)
    buf_y1 = y+dt*buf_f0
    CALL ad(t+dt,ystar,buf_y1,buf_f1)
    res=y+0.5*(buf_f0+buf_f1)*dt
    t=t+dt
  END SUBROUTINE ad_step

  !> Routine to perform a simultaneous integration step (RK4 algorithm) of the nonlinear and adjoint together. The incremented time is returned.
  !> @param y Model variable at time t
  !> @param deltay Perturbation at time t
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param ynew Model variable at time t+dt
  !> @param deltaynew Perturbation at time t+dt
  SUBROUTINE evolve_ad_step(y,deltay,t,dt,ynew,deltaynew)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y,deltay
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: ynew,deltaynew

    CALL tendencies(t,y,buf_f0)
    CALL ad(t,y,deltay,buf_f00)

    buf_y1 = y + dt*buf_f0
    buf_y11 = deltay + dt*buf_f00

    CALL tendencies(t+dt,buf_y1,buf_f1)
    CALL ad(t+dt,buf_y1,buf_y11,buf_f11)
    
    t=t+dt
    ynew=y+0.5*(buf_f0+buf_f1)*dt
    deltaynew=deltay+0.5*(buf_f00+buf_f11)*dt
  END SUBROUTINE evolve_ad_step


  !-----------------------------------------------------!
  !                                                     !
  ! Tangent linear model integrator                     !
  !                                                     !
  !-----------------------------------------------------!

  !> Routine to perform an integration step (Heun algorithm) of the tangent linear model. The incremented time is returned.
  !> @param y Initial point.
  !> @param ystar Adjoint model at the point ystar.
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param res Final point after the step.
  SUBROUTINE tl_step(y,ystar,t,dt,res)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y,ystar
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res

    CALL tl(t,ystar,y,buf_f0)
    buf_y1 = y+dt*buf_f0
    CALL tl(t+dt,ystar,buf_y1,buf_f1)
    res=y+0.5*(buf_f0+buf_f1)*dt
    t=t+dt
  END SUBROUTINE tl_step

  !> Routine to perform a simultaneous integration step (RK4 algorithm) of the nonlinear and the tangent linear model together. The incremented time is returned.
  !> @param y Model variable at time t
  !> @param deltay Perturbation at time t
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param ynew Model variable at time t+dt
  !> @param deltaynew Perturbation at time t+dt
  SUBROUTINE evolve_tl_step(y,deltay,t,dt,ynew,deltaynew)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y,deltay
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: ynew,deltaynew

    CALL tendencies(t,y,buf_f0)
    CALL ad(t,y,deltay,buf_f00)

    buf_y1 = y + dt*buf_f0
    buf_y11 = deltay + dt*buf_f00

    CALL tendencies(t+dt,buf_y1,buf_f1)
    CALL ad(t+dt,buf_y1,buf_y11,buf_f11)
    
    t=t+dt
    ynew=y+0.5*(buf_f0+buf_f1)*dt
    deltaynew=deltay+0.5*(buf_f00+buf_f11)*dt
  END SUBROUTINE evolve_tl_step


  !> Routine to perform a simultaneously an integration step (Heun algorithm) of the nonlinear and computes the Heun tangent linear propagator. The boolean variable adjoint allows for an adjoint forward integration. The incremented time is returned.
  !> @param y Model variable at time t
  !> @param propagator Propagator at time t
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param ynew Model variable at time t+dt
  !> @param adjoint If true, compute the propagator of the adjoint model (AD) instead of the tangent one (TL)
  SUBROUTINE prop_step(y,propagator,t,dt,ynew,adjoint)
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    LOGICAL, INTENT(IN) :: adjoint
    REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(OUT) :: propagator
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: ynew
 
    CALL tendencies(t,y,buf_f0)
    buf_j1=jacobian_mat(y)

    buf_y1 = y + dt*buf_f0

    CALL tendencies(t+dt,buf_y1,buf_f1)
    buf_j2=jacobian_mat(buf_y1)
    
    buf_j1h=buf_j1
    buf_j2h=buf_j2
    CALL dgemm ('n', 'n', ndim, ndim, ndim, dt, buf_j2, ndim,buf_j1h, ndim,1.0d0, buf_j2h, ndim)
     
    ynew=y  + dt/2.0d0*(buf_f0 + buf_f1)
    IF (adjoint) THEN
            propagator=one - dt/2.0d0*(buf_j1h + buf_j2h)
    ELSE
            propagator=one + dt/2.0d0*(buf_j1h + buf_j2h)
    END IF        
    t=t+dt
   
  END SUBROUTINE prop_step


  
END MODULE tl_ad_integrator
