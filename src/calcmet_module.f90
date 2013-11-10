module calcmet_module
  ! 
  ! DESCRIPTION
  !  this module contains the functions to calculate meteorological variables
  !
  use type_module
  use constant_module, only : grav=>earth_gravity, air_cp, gascon=>air_rd, isa_gamma
  use p2p_module, only : extrapolate_Ts
  use ip_module, only : searchidx
  implicit none

  private
  public :: calcmet_q, calcmet_rh, calcmet_ps, calcmet_z, calcmet_rh_fromTd, calcmet_msl

  interface calcmet_q
     module procedure calcmet_q_0d, calcmet_q_1d0d, calcmet_q_2d0d, &
          &           calcmet_q_3d1d, calcmet_q_3d3d
  end interface calcmet_q
  interface calcmet_rh
     module procedure calcmet_rh_0d, calcmet_rh_1d0d, calcmet_rh_2d0d, &
          &           calcmet_rh_3d1d, calcmet_rh_3d3d
  end interface calcmet_rh
  interface calcmet_z
     module procedure calcmet_z_3d1d,calcmet_z_3d3d
  end interface calcmet_z
  interface calcmet_rh_fromTd
     module procedure calcmet_rh_fromTd_0d, calcmet_rh_fromTd_3d
  end interface calcmet_rh_fromTd
  interface calcmet_msl
     module procedure calcmet_msl_1d1d, calcmet_msl_3d3d
  end interface calcmet_msl
  
contains
  !!
  !
  ! FUNCTION calcmet_q(T, rh, plev) result (q)
  !
  ! DESCRIPTION
  !  compulte specific humidity using relative humidity
  ! 
  ! ARGUMENTS
  !  INPUT:
  !    T    : temperature       [K]
  !    rh   : relative humdity  [%]
  !    plev : pressure          [Pa]
  !  RETURN:
  !    q    : sprecific humidity [Kg/Kg]
  !
  ! NOTE
  !  if T is array, the same dimension array 'q' is retured  
  !
  !!
  function calcmet_q_0d(T, rh, plev) result (q)    
    real(kind=sp), intent(in) :: T, rh
    real(kind=sp), intent(in) :: plev       !Pa
    real(kind=sp) :: q

    real(kind=sp), parameter :: eps = 0.622 ! Rd/Rv
    real(kind=sp) :: es, e

    es = exp(19.482 - 4303.4 / (T-29.65)) ! in hPa JMA/WMO
    e = es*rh ! x 100 (hPa->Pa) / 100. (%->no dimension)
    q = eps*e / (plev - (1.0-eps)*e)

  end function calcmet_q_0d
  function calcmet_q_1d0d(T, rh, plev) result (q)
    real(kind=sp), dimension(:), intent(in) :: T, rh
    real(kind=sp), intent(in) :: plev
    real(kind=sp), dimension(size(T,1)) :: q

    integer(kind=i4b) :: i

    do i=1, size(T,1)
       q(i) = calcmet_q_0d(T(i), rh(i), plev)
    end do

  end function calcmet_q_1d0d
  function calcmet_q_2d0d(T, rh, plev) result (q)
    real(kind=sp), dimension(:,:), intent(in) :: T, rh
    real(kind=sp), intent(in) :: plev
    real(kind=sp), dimension(size(T,1), size(T,2)) :: q

    integer(kind=i4b) :: i, j

    do j=1, size(T,2)
      do i=1, size(T,1)
         q(i,j) = calcmet_q_0d(T(i,j), rh(i,j), plev)
      end do
    end do

  end function calcmet_q_2d0d
  function calcmet_q_3d1d(T, rh, plev) result (q)
    real(kind=sp), dimension(:,:,:), intent(in) :: T, rh
    real(kind=sp), dimension(:), intent(in) :: plev
    real(kind=sp), dimension(size(T,1), size(T,2), size(T,3)) :: q

    integer(kind=i4b) :: i, j, k

    do k = 1, size(T,3)
       do j=1, size(T,2)
          do i=1, size(T,1)
             q(i,j,k) = calcmet_q_0d(T(i,j,k), rh(i,j,k), plev(k))
          end do
       end do
    end do

  end function calcmet_q_3d1d
  function calcmet_q_3d3d(T, rh, plev) result (q)
    real(kind=sp), dimension(:,:,:), intent(in) :: T, rh, plev
    real(kind=sp), dimension(size(T,1), size(T,2), size(T,3)) :: q

    integer(kind=i4b) :: i, j, k

    do k = 1, size(T,3)
       do j=1, size(T,2)
          do i=1, size(T,1)
             q(i,j,k) = calcmet_q_0d(T(i,j,k), rh(i,j,k), plev(i,j,k))
             if( q(i,j,k)<0 ) then
             end if
          end do
       end do
    end do

  end function calcmet_q_3d3d

  !!
  ! 
  ! FUNCTION calcmet_rh(T, q, plev) result (rh)
  !
  ! DESCRIPTION
  !  compute relative humidity using specific humidity
  !
  ! ARGUMENTS
  !  INPUT:
  !    T    : temperature      [K]
  !    q    : specfic humidity [%]
  !    plev : pressure         [Pa]
  !  OUTPUT:
  !    rh   : relative humidity [%]
  !
  ! NOTE
  !  if T is array, the same dimension array 'rh' is retured  
  !
  !!
  function calcmet_rh_0d(T, q, plev) result (rh)
    real(kind=sp), intent(in) :: T, q
    real(kind=sp), intent(in) :: plev
    real(kind=sp) :: rh

    real(kind=sp), parameter :: eps = 0.622 ! Rd/Rv
    real(kind=sp) :: es, e

    es = exp(19.482 - 4303.4 / (T-29.65)) ! in hPa JMA/WMO
    e = q*plev / (eps + (1-eps)*q)
    rh = e/es ! x 100 (hPa->Pa) / 100. (%->no dimension)

  end function calcmet_rh_0d
  function calcmet_rh_1d0d(T, q, plev) result (rh)
    real(kind=sp), dimension(:), intent(in) :: T, q
    real(kind=sp), intent(in) :: plev
    real(kind=sp), dimension(size(T,1)) :: rh

    integer(kind=i4b) :: i

    do i = 1, size(T,1)
       rh(i) = calcmet_rh_0d(T(i), q(i), plev)
    end do

  end function calcmet_rh_1d0d
  function calcmet_rh_2d0d(T, q, plev) result (rh)
    real(kind=sp), dimension(:,:), intent(in) :: T, q
    real(kind=sp), intent(in) :: plev
    real(kind=sp), dimension(size(T,1), size(T,2)) :: rh

    integer(kind=i4b) :: i, j

    do j = 1, size(T,2)
       do i = 1, size(T,1)
          rh(i,j) = calcmet_rh_0d(T(i,j), q(i,j), plev)
       end do
    end do

  end function calcmet_rh_2d0d
  function calcmet_rh_3d1d(T, q, plev) result (rh)
    real(kind=sp), dimension(:,:,:), intent(in) :: T, q
    real(kind=sp), dimension(:), intent(in) :: plev
    real(kind=sp), dimension(size(T,1), size(T,2), size(T,3)) :: rh
    integer(kind=i4b) :: i, j, k

    do k = 1, size(T,3)
       do j = 1, size(T,2)
          do i = 1, size(T,1)
             rh(i,j,k) = calcmet_rh_0d(T(i,j,k), q(i,j,k), plev(k))
          end do
       end do
    end do

  end function calcmet_rh_3d1d
  function calcmet_rh_3d3d(T, q, plev) result (rh)
    real(kind=sp), dimension(:,:,:), intent(in) :: T, q
    real(kind=sp), dimension(:,:,:), intent(in) :: plev
    real(kind=sp), dimension(size(T,1), size(T,2), size(T,3)) :: rh
    integer(kind=i4b) :: i, j, k

    do k = 1, size(T,3)
       do j = 1, size(T,2)
          do i = 1, size(T,1)
             rh(i,j,k) = calcmet_rh_0d(T(i,j,k), q(i,j,k), plev(i,j,k))
          end do
       end do
    end do

  end function calcmet_rh_3d3d

  !!
  ! FUNCTION calcmet_ps(p, T, z, zs) result(ps)
  !
  ! DESCRIPTION
  !  compute surface pressure using hydrostatic equation
  !
  ! ARGUMENTS
  !  INPUT:
  !    p  : pressure             [Pa]
  !    T  : temperature          [K]
  !    z  : geopotential height  [m]
  !    zs : geometric height     [m]
  !  RETURN:
  !    ps : surface pressure     [Pa]
  !
  !!
  function calcmet_ps(p, T, z, zs) result(ps)
    real(kind=sp), dimension(:,:,:), intent(in) :: p, T, z
    real(kind=sp), dimension(:,:), intent(in) :: zs
    real(kind=sp), dimension(size(zs,1), size(zs,2)) :: ps

    real(kind=sp), parameter :: ggg = grav/gascon/isa_gamma
    integer(kind=i4b), dimension(1) :: kk
    integer(kind=i4b) :: i, j, k

    do j=1, size(zs,2)
      do i=1, size(zs, 1)
        kk = minloc(abs(zs(i,j)-z(i,j,:)),mask=zs(i,j)<z(i,j,:))
        k = kk(1)
        ps(i,j) = p(i,j,k) * (1.0 - isa_gamma*(zs(i,j)-z(i,j,k))/T(i,j,k))**ggg
      end do
    end do
  end function calcmet_ps

  !!
  ! FUNCTION calcmet_z(T, zs, ps, pl) result(z)
  !
  ! DESCRIPTION
  !  compute geopotential height using hydrostatic equation
  !
  ! Argument
  !  INPUT:
  !    T  : temperature      [K]
  !    zs : geometirc height [m]
  !    ps : surface pressure [Pa]
  !    pl : pressure         [Pa]
  !
  !  Return
  !    z  : geopotential height [m]
  !
  !!
  function calcmet_z_3d1d(T,zs,ps,pl) result(z)
    real(kind=sp), dimension(:,:,:), intent(in) :: T 
    real(kind=sp), dimension(:,:), intent(in) :: zs,ps 
    real(kind=sp), dimension(:), intent(in) :: pl    

    real(kind=sp), parameter :: gg = gascon/grav
    real(kind=sp), dimension(size(T,1),size(T,2),size(T,3)) :: z

    integer(kind=i4b) :: k, n

    n = size(pl)
    z(:,:,1) = zs(:,:) - gg*T(:,:,1)*log(pl(1)/ps(:,:))

    do k=2, n
      z(:,:,k) = z(:,:,k-1) + 0.5*gg*(T(:,:,k-1)+T(:,:,k))*(log(pl(k-1))-log(pl(k)))
    end do

  end function calcmet_z_3d1d
  function calcmet_z_3d3d(T,zs,ps,pl) result(z)
    real(kind=sp), dimension(:,:,:), intent(in) :: T 
    real(kind=sp), dimension(:,:), intent(in) :: zs,ps 
    real(kind=sp), dimension(:,:,:), intent(in) :: pl 

    real(kind=sp), parameter :: gg = gascon/grav
    real(kind=sp), dimension(size(T,1),size(T,2),size(T,3)) :: z

    integer(kind=i4b) :: k, n

    n = size(pl,3)
    z(:,:,1) = zs(:,:) - gg*T(:,:,1)*log(pl(:,:,1)/ps(:,:))
    do k=2, n
      z(:,:,k) = z(:,:,k-1) + 0.5*gg*(T(:,:,k-1)+T(:,:,k))*(log(pl(:,:,k-1))-log(pl(:,:,k)))
    end do

  end function calcmet_z_3d3d

  !!
  ! FUNCTION calcmet_rh_fromTd(T, Td) result(rh)
  !
  ! DESCRIPTION
  !  compute relative humidity using dew point temperature
  !
  ! ARGUMENTS
  !  INPUT:
  !    T  : temperature [K]
  !    Td : dew pointe temperature [K]
  !  RETURN
  !    rh : relative humidity [%]
  !
  !!
  function calcmet_rh_fromTd_0d(T, Td) result (rh)
    real(kind=sp), intent(in) :: T, Td
    real(kind=sp) :: rh

    real(kind=sp), parameter :: eps = 0.622 ! Rd/Rv
    real(kind=sp) :: es, e

    es = exp(19.482 - 4303.4 / (T-29.65))  ! in hPa JMA/WMO
    e  = exp(19.482 - 4303.4 / (Td-29.65)) ! in hPa JMA/WMO
    rh = e/es * 100 

  end function calcmet_rh_fromTd_0d
  function calcmet_rh_fromTd_3d(T, Td) result (rh)
    real(kind=sp), dimension(:,:,:), intent(in) :: T, Td
    real(kind=sp), dimension(size(T,1), size(T,2), size(T,3)) :: rh
    integer(kind=i4b) :: i, j, k

    do k = 1, size(T,3)
       do j = 1, size(T,2)
          do i = 1, size(T,1)
             rh(i,j,k) = calcmet_rh_fromTd_0d(T(i,j,k), Td(i,j,k))
          end do
       end do
    end do

  end function calcmet_rh_fromTd_3d

  !!
  ! FUNCTION calcmet_msl() result (msl)
  !
  ! DESCRIPTION
  !  compute mean sea level pressure using hydrostatic equation
  !
  ! ARGUMENTS
  !  INPUT:
  !    T    : temperature             [K]
  !    zs   : geopmetoric height      [m]
  !    ps   : surface pressure        [Pa]
  !    plev : pressure level          [Pa]
  !  OUTPUT:
  !    msl  : mean sea level pressure [Pa]
  !
  !!
  function calcmet_msl_1d1d(T,zs,ps,plev) result (msl)
    real(kind=sp), intent(in) :: T(:),plev(:)
    real(kind=sp), intent(in) :: zs,ps
    real(kind=sp), parameter  :: dTdz=isa_gamma
    real(kind=sp)             :: msl
    integer(kind=i4b) :: kl
    real(kind=sp)     :: Ts

    kl = searchidx(plev,ps,-1)                     !plev(kl-1)>ps>plev(kl)
    Ts = extrapolate_Ts(T(kl),plev(kl),ps)
    msl = ps*(1.+dTdz*zs/Ts )**(grav/gascon/dTdz)

  end function calcmet_msl_1d1d
  function calcmet_msl_3d3d(T,zs,ps,plev) result (msl)
    real(kind=sp), intent(in) :: T(:,:,:),plev(:,:,:)
    real(kind=sp), intent(in) :: zs(:,:),ps(:,:)
    real(kind=sp) :: msl(size(zs,1) ,size(zs,2))
    integer(kind=i4b) :: i,j

    do j = 1, size(zs,2)
       do i = 1, size(zs,1)
          msl(i,j) = calcmet_msl_1d1d(T(i,j,:),zs(i,j),ps(i,j),plev(i,j,:))
       end do
    end do

  end function calcmet_msl_3d3d

end module calcmet_module
