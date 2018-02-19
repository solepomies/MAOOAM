
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
!>  This module actually contains the RK4 algorithm routines.
!>  The user can modify it according to its preferred integration scheme.
!>  For higher-order schemes, additional bufers will probably have to be defined.
!                                                                           
!---------------------------------------------------------------------------



MODULE tl_ad_integrator

  USE util, only: init_one
  USE params, only: ndim
  USE tensor, only: sparse_mul3
  USE aotensor_def, only: aotensor

  USE tl_ad_tensor , only: ad,tl,jacobian_mat
  IMPLICIT NONE

  PRIVATE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y1 !< Buffer to hold the intermediate position of the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_kA !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_kB !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y11 !< Buffer to hold the intermediate position of the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_kC !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_kD !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j1 !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j2 !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j3 !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j4 !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j1h !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j2h !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j3h !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: buf_j4h !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_kAA !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_kBB !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: one        !< unit matrix 

    
  PUBLIC :: ad_step, init_tl_ad_integrator, tl_step, evolve_ad_step, evolve_tl_step, prop_step

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

  !> Routine to initialise the TL-AD integration bufers.
  SUBROUTINE init_tl_ad_integrator
    INTEGER :: AllocStat
    ALLOCATE(buf_j1h(ndim,ndim),buf_j2h(ndim,ndim),buf_j3h(ndim,ndim),&
         &buf_j4h(ndim,ndim),buf_j1(ndim,ndim),buf_j2(ndim,ndim),&
         &buf_j3(ndim,ndim),buf_j4(ndim,ndim),one(ndim,ndim),buf_y11(0:ndim),&
         &buf_y1(0:ndim),buf_kA(0:ndim),buf_kB(0:ndim),buf_kC(0:ndim),buf_kD(0:ndim), &
         &buf_kAA(0:ndim),buf_kBB(0:ndim) ,STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    CALL init_one(one)
  END SUBROUTINE init_tl_ad_integrator


  !-----------------------------------------------------!
  !                                                     !
  ! Adjoint model integrator                            !
  !                                                     !
  !-----------------------------------------------------!


  

  !> Routine to perform an integration step (RK4 algorithm) of the adjoint model. The incremented time is returned.
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

    CALL ad(t,ystar,y,buf_kA)
    buf_y1 = y+0.5*dt*buf_kA
    CALL ad(t+0.5*dt,ystar,buf_y1,buf_kB)
    buf_y1 = y+0.5*dt*buf_kB
    buf_kA = buf_kA+2*buf_kB
    CALL ad(t+0.5*dt,ystar,buf_y1,buf_kB)
    buf_y1 = y+0.5*dt*buf_kB
    buf_kA = buf_kA+2*buf_kB
    CALL ad(t+dt,ystar,buf_y1,buf_kB)
    buf_kA = buf_kA+buf_kB
    res=y+buf_kA*dt/6
    t=t+dt
  END SUBROUTINE ad_step

  !> Routine to perform a simultaneous integration step (RK4 algorithm) of the nonlinear and adjoint model together. The incremented time is returned.
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

    CALL tendencies(t,y,buf_kA)
    CALL ad(t,y,deltay,buf_kAA)

    buf_y1 = y + 0.5*dt*buf_kA
    buf_y11 = deltay + 0.5*dt*buf_kAA

    CALL tendencies(t+0.5*dt,buf_y1,buf_kB)
    CALL ad(t+0.5*dt,buf_y1,buf_y11,buf_kBB)
    
    buf_y1 = y + 0.5*dt*buf_kB
    buf_y11 = deltay + 0.5*dt*buf_kBB

    buf_kA = buf_kA + 2*buf_kB
    buf_kAA = buf_kAA + 2*buf_kBB
    
    CALL tendencies(t+0.5*dt,buf_y1,buf_kB)
    CALL ad(t+0.5*dt,buf_y1,buf_y11,buf_kBB)
    
    buf_y1 = y + dt*buf_kB
    buf_y11 = deltay + dt*buf_kBB
    
    buf_kA = buf_kA + 2*buf_kB
    buf_kAA = buf_kAA + 2*buf_kBB
    
    CALL tendencies(t+dt,buf_y1,buf_kB)
    CALL ad(t+dt,buf_y1,buf_y11,buf_kBB)

    buf_kA = buf_kA + buf_kB
    buf_kAA = buf_kAA + buf_kBB
    
    t=t+dt
    ynew=y+buf_kA*dt/6
    deltaynew=deltay+buf_kAA*dt/6
  END SUBROUTINE evolve_ad_step


  !-----------------------------------------------------!
  !                                                     !
  ! Tangent linear model integrator                     !
  !                                                     !
  !-----------------------------------------------------!

  !> Routine to perform an integration step (RK4 algorithm) of the tangent linear model. The incremented time is returned.
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

    CALL tl(t,ystar,y,buf_kA)
    buf_y1 = y+0.5*dt*buf_kA
    CALL tl(t+0.5*dt,ystar,buf_y1,buf_kB)
    buf_y1 = y+0.5*dt*buf_kB
    buf_kA = buf_kA+2*buf_kB
    CALL tl(t+0.5*dt,ystar,buf_y1,buf_kB)
    buf_y1 = y+0.5*dt*buf_kB
    buf_kA = buf_kA+2*buf_kB
    CALL tl(t+dt,ystar,buf_y1,buf_kB)
    buf_kA = buf_kA+buf_kB
    res=y+buf_kA*dt/6
    t=t+dt
  END SUBROUTINE tl_step

  !> Routine to perform a simultaneous integration step (RK4 algorithm) of the nonlinear and tangent linear model togheter. The incremented time is returned.
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

    CALL tendencies(t,y,buf_kA)
    CALL tl(t,y,deltay,buf_kAA)

    buf_y1 = y + 0.5*dt*buf_kA
    buf_y11 = deltay + 0.5*dt*buf_kAA

    CALL tendencies(t+0.5*dt,buf_y1,buf_kB)
    CALL tl(t+0.5*dt,buf_y1,buf_y11,buf_kBB)
    
    buf_y1 = y + 0.5*dt*buf_kB
    buf_y11 = deltay + 0.5*dt*buf_kBB

    buf_kA = buf_kA + 2*buf_kB
    buf_kAA = buf_kAA + 2*buf_kBB
    
    CALL tendencies(t+0.5*dt,buf_y1,buf_kB)
    CALL tl(t+0.5*dt,buf_y1,buf_y11,buf_kBB)
    
    buf_y1 = y + dt*buf_kB
    buf_y11 = deltay + dt*buf_kBB
    
    buf_kA = buf_kA + 2*buf_kB
    buf_kAA = buf_kAA + 2*buf_kBB
    
    CALL tendencies(t+dt,buf_y1,buf_kB)
    CALL tl(t+dt,buf_y1,buf_y11,buf_kBB)

    buf_kA = buf_kA + buf_kB
    buf_kAA = buf_kAA + buf_kBB
    
    t=t+dt
    ynew=y+buf_kA*dt/6
    deltaynew=deltay+buf_kAA*dt/6
  END SUBROUTINE evolve_tl_step

  !> Routine to perform a simultaneously an integration step (RK4 algorithm) of the nonlinear and computes the RK4 tangent linear propagator. The boolean variable adjoint allows for an adjoint forward integration. The incremented time is returned.
  !> @param y Model variable at time t
  !> @param propagator Propagator at time t
  !> @param t Actual integration time
  !> @param dt Integration timestep
  !> @param ynew Model variable at time t+dt
  !> @param adjoint If true, compute the propagator of the adjoint model (AD) instead of the tangent one (TL)
  SUBROUTINE prop_step(y,propagator,t,dt,ynew,adjoint)
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    LOGICAL, INTENT(IN) :: adjoint
    REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(OUT) :: propagator
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: ynew
 
    CALL tendencies(t,y,buf_kA)
    buf_j1=jacobian_mat(y)

    buf_y1 = y + 0.5*dt*buf_kA

    CALL tendencies(t+0.5*dt,buf_y1,buf_kB)
    buf_j2=jacobian_mat(buf_y1)
    
    buf_y1 = y + 0.5*dt*buf_kB

    CALL tendencies(t+0.5*dt,buf_y1,buf_kC)
    buf_j3=jacobian_mat(buf_y1)

    buf_y1 = y + dt*buf_kC
    CALL tendencies(t+dt,buf_y1,buf_kD)
    buf_j4=jacobian_mat(buf_y1)
    
    buf_j1h=buf_j1
    buf_j2h=buf_j2
    buf_j3h=buf_j3
    buf_j4h=buf_j4
    call dgemm ('n', 'n', ndim, ndim, ndim, dt/2.0d0, buf_j2, ndim,buf_j1h, ndim,1.0d0, buf_j2h, ndim)
    call dgemm ('n', 'n', ndim, ndim, ndim, dt/2.0d0, buf_j3, ndim,buf_j2h, ndim,1.0d0, buf_j3h, ndim)
    call dgemm ('n', 'n', ndim, ndim, ndim, dt      , buf_j4, ndim,buf_j3h, ndim,1.0d0, buf_j4h, ndim)
     
    ynew=y  + dt/6.0d0*(buf_kA + 2.0d0*buf_kB + 2.0d0*buf_kC + buf_kD)
    IF (adjoint) THEN
            propagator=one - dt/6.0d0*(buf_j1h + 2.0d0*buf_j2h + 2.0d0*buf_j3h + buf_j4h)
    ELSE
            propagator=one + dt/6.0d0*(buf_j1h + 2.0d0*buf_j2h + 2.0d0*buf_j3h + buf_j4h)
    END IF        
    t=t+dt
   
  END SUBROUTINE prop_step

END MODULE tl_ad_integrator
