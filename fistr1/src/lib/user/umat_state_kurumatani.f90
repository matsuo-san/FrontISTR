!-------------------------------------------------------------------------------
! Copyright (c) 2019 FrontISTR Commons
! This software is released under the MIT License, see LICENSE.txt
!-------------------------------------------------------------------------------
!>  \brief   This subroutine read in used-defined material properties
!>  tangent
module mUmat
  implicit none

  integer, parameter, private :: kreal = kind(0.0d0)

contains

  subroutine mate( matl, strain, D, fstat )

    implicit none

    real(kind=kreal), intent(in)  :: strain(6) !< Green-Lagrangen strain
    real(kind=kreal), intent(in)  :: matl(:)   !< material properties
    real(kind=kreal), intent(out) :: D(:,:)    !< strain-stress relation
    real(kind=kreal), intent(inout)  :: fstat(:) 
    real(kind=kreal) :: EE, PP

    real(kind=kreal) :: damage
    real(kind=kreal) :: sdev(6), I1, J2
    real(kind=kreal) :: epse, eps0, gf
    real(kind=kreal) :: kcomp

    D(:,:) = 0.d0
    damage = 0.d0

    EE = matl(1)
    PP = matl(2)

    !
    kcomp = matl(5)
    !
    I1      = strain(1)+strain(2)+strain(3)
    sdev(:) = strain(:)-1.d0/3.d0*I1
    J2 = 0.5d0*( sdev(1)**2.d0+sdev(2)**2.d0+sdev(3)**2.d0 &
              &  + 2.d0*(sdev(4)**2.d0+sdev(5)**2.d0+sdev(6)**2.d0) )
    epse = (kcomp-1.d0)/(2.d0*kcomp*(1.d0-2.d0*PP))
    epse = epse &
         &  + 1.d0/(2.d0*kcomp)*sqrt(((kcomp-1)*I1/(1.d0-2.d0*PP))**2.d0+12.d0*kcomp*J2/(1.d0+PP)**2.d0)
    eps0 = matl(3)
    !
    epse = max( epse, fstat(1) )
    fstat(1) = epse
    !
    !epse = strain(3)
    !
    if( epse>eps0 ) then
      gf     = matl(4)
      damage = 1.d0-(eps0/epse)*exp(-1.d0*EE*eps0*(epse-eps0)/gf)
      damage = min( damage,0.9995d0 )
      !write(*,*) damage, epse, strain(:)
      !write(*,*) damage, exp(-1.d0*EE*eps0*(epse-eps0)/gf), (-1.d0*EE*eps0*(epse-eps0)/gf), epse, I1, J2!, ( (kcomp-1.d0)/( 2.0d0 * ( 1.0d0 - 2.0d0 * PP ) ) )
      !write(*,*) damage, exp(-1.d0*EE*eps0*(epse-eps0)/gf), (-1.d0*EE*eps0*(epse-eps0)/gf), epse, I1, J2, (kcomp-1.d0)/( 2.d0 * kcomp * (1.d0-2.d0*PP) )
      !1.d0/(2.d0*kcomp)*sqrt(((kcomp-1)*I1/(1.d0-2.d0*PP))**2.d0+12.d0*kcomp*J2/(1.d0+PP)**2.d0)
      !damage = min(damage,0.995d0)
      !damage = min(damage,0.99d0)
      !EE = EE *(1.d0-damage)
    ! damage = max(damage,1.d-3)
    ! damage = 1.d0-damage
    endif
    !
    D(1,1)=EE*(1.d0-PP)/(1.d0-2.d0*PP)/(1.d0+PP)
    D(1,2)=EE*PP/(1.d0-2.d0*PP)/(1.d0+PP)
    D(1,3)=D(1,2)
    D(2,1)=D(1,2)
    D(2,2)=D(1,1)
    D(2,3)=D(1,2)
    D(3,1)=D(1,3)
    D(3,2)=D(2,3)
    D(3,3)=D(1,1)
    D(4,4)=EE/(1.d0+PP)*0.5d0
    D(5,5)=EE/(1.d0+PP)*0.5d0
    D(6,6)=EE/(1.d0+PP)*0.5d0

    D(:,:) = (1.d0-damage)*D(:,:)

  end subroutine mate

  !> This subroutine calculates constitutive matrix
  subroutine uMatlMatrix( mname, matl, strain, stress, fstat, D  &
      , dtime, ttime, temperature )
    character(len=*), intent(in)  :: mname     !< material name
    real(kind=kreal), intent(in)  :: matl(:)   !< material properties
    real(kind=kreal), intent(in)  :: strain(6) !< Green-Lagrangen strain
    real(kind=kreal), intent(in)  :: stress(6) !< 2nd Piola-Kirchhiff stress tensor
    real(kind=kreal), intent(inout)  :: fstat(:)  !< state variables
    real(kind=kreal), intent(out) :: D(:,:)    !< strain-stress relation
    real(kind=kreal), intent(in)  :: dtime     !< time increment
    real(kind=kreal), intent(in)  :: ttime     !< total time at the start of the current increment
    real(kind=kreal), optional    :: temperature !< temprature

    !
    call mate( matl(:), strain, D, fstat )

  end subroutine

  !> This subroutine calculate strain and stress increment
  subroutine uUpdate(  mname, matl, strain, stress, fstat, dtime, ttime, temperature )
    character(len=*), intent(in)    :: mname      !< material name
    real(kind=kreal), intent(in)    :: matl(:)    !< material properties
    real(kind=kreal), intent(in)    :: strain(6)  !< strain
    real(kind=kreal), intent(inout) :: stress(6)  !< 2nd Piola-Kirchhiff stress tensor
    real(kind=kreal), intent(inout) :: fstat(:)   !< state variables
    real(kind=kreal), intent(in)    :: dtime     !< time increment
    real(kind=kreal), intent(in)    :: ttime     !< total time at the start of the current increment
    real(kind=kreal), optional      :: temperature !< temperature

    real(kind=kreal) :: D(6,6)

    !call uElasticMatrix( matl(:2), strain, D )
    call mate( matl(:), strain, D, fstat )

    if( fstat(1) < matl(3) ) then
      stress = matmul( D, strain ) + stress
    else
      !stress = -1.d0*matmul( D, strain ) + stress
      stress = matmul( D, strain )
    endif

  end subroutine

end module
