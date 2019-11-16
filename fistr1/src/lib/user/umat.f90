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

  subroutine calc_equivalent_strain( matl, strain, epsmax, epse )

    implicit none

    real(kind=kreal), intent(in)  :: strain(6) !< Green-Lagrangen strain
    real(kind=kreal), intent(in)  :: matl(:)   !< material properties
    real(kind=kreal), intent(in)  :: epsmax
    real(kind=kreal), intent(out) :: epse
    !
    real(kind=kreal) :: sdev(6), I1, J2
    real(kind=kreal) :: PP
    real(kind=kreal) :: kcomp

    PP    = matl(2)
    kcomp = matl(5)
    !
    I1      = strain(1)+strain(2)+strain(3)
    sdev(:) = strain(:)-1.d0/3.d0*I1
    J2 = 0.5d0*( sdev(1)**2.d0+sdev(2)**2.d0+sdev(3)**2.d0 &
              &  + 2.d0*(sdev(4)**2.d0+sdev(5)**2.d0+sdev(6)**2.d0) )
    !
    epse = (kcomp-1.d0)/(2.d0*kcomp*(1.d0-2.d0*PP))
    epse = epse &
         &  + 1.d0/(2.d0*kcomp)*sqrt(((kcomp-1)*I1/(1.d0-2.d0*PP))**2.d0+12.d0*kcomp*J2/(1.d0+PP)**2.d0)
    !
    epse = max( epse, epsmax )

  endsubroutine calc_equivalent_strain
  !
  subroutine calc_damage( matl,epse,length,damage )

    implicit none

    real(kind=kreal), intent(in)  :: matl(:)   !< material properties
    real(kind=kreal), intent(in)  :: epse
    real(kind=kreal), intent(in)  :: length
    real(kind=kreal), intent(out) :: damage
    real(kind=kreal) :: eps0, gf, EE

    eps0 = matl(3)

    if( epse>eps0 ) then
      EE     = matl(1)
      gf     = matl(4)
      damage = 1.d0-(eps0/epse)*exp(-1.d0*EE*length*eps0*(epse-eps0)/gf)
      damage = min( damage,0.9999d0 )
    else
      damage = 0.d0
    endif
    !
  endsubroutine calc_damage
  !
  subroutine calc_elastic_matrix( matl, D )

    implicit none

    real(kind=kreal), intent(in)  :: matl(:)   !< material properties
    real(kind=kreal), intent(out) :: D(:,:)    !< strain-stress relation
    real(kind=kreal) :: EE, PP

    D(:,:) = 0.d0
    !
    EE = matl(1)
    PP = matl(2)
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

  endsubroutine calc_elastic_matrix
  !
  !> This subroutine calculates constitutive matrix
  subroutine uMatlMatrix( mname, matl, strain, stress, fstat, D  &
      , dtime, ttime, temperature )
    character(len=*), intent(in)  :: mname     !< material name
    real(kind=kreal), intent(in)  :: matl(:)   !< material properties
    real(kind=kreal), intent(in)  :: strain(6) !< Green-Lagrangen strain
    real(kind=kreal), intent(in)  :: stress(6) !< 2nd Piola-Kirchhiff stress tensor
    real(kind=kreal), intent(in)  :: fstat(:)  !< state variables
    real(kind=kreal), intent(out) :: D(:,:)    !< strain-stress relation
    real(kind=kreal), intent(in)  :: dtime     !< time increment
    real(kind=kreal), intent(in)  :: ttime     !< total time at the start of the current increment
    real(kind=kreal), optional    :: temperature !< temprature

    real(kind=kreal) :: epse
    real(kind=kreal) :: damage

    call calc_elastic_matrix( matl, D )
    call calc_equivalent_strain( matl, strain, fstat(13), epse )
    call calc_damage( matl,epse,fstat(16),damage )

    D(:,:) = (1.d0-damage)*D(:,:)

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

    real(kind=kreal) :: epse
    real(kind=kreal) :: damage

    !
    ! fstat(1:6)  :: strain at current time (trial)
    ! fstat(7:12) :: strain at previous time
    ! fstat(13)   :: maximum equiv strain at current time (trial)
    ! fstat(14)   :: maximum equiv strain at previous time
    ! fstat(15)   :: previous time
    !

    !
    ! update when time incremented
    if( ttime-fstat(15) > 1e-8 ) then
      fstat(7:12) = fstat(1:6)
      fstat(14)  = fstat(13)
      fstat(15)   = ttime
    !
    else
      fstat(1:6) = fstat(7:12)
      fstat(13)  = fstat(14)
    endif

    fstat(1:6) = fstat(1:6) + strain(:)
    !
    call calc_elastic_matrix( matl, D )
    call calc_equivalent_strain( matl, fstat(1:6), fstat(13), epse )
    call calc_damage( matl,epse,fstat(16),damage )
    !
    D(:,:) = (1.d0-damage)*D(:,:)
    !
    stress    = matmul( D, fstat(1:6) )
    fstat(13) = epse
    !
  end subroutine

end module
