module gfs_module
  !
  ! DESCRIPTION
  !  this module contains the wrapper routines to read & write gfs sigma file
  !
  use type_module
  use constant_module, only : PI=>math_pi  
  use sigio_module
  use sigio_r_module
  implicit none
  
  private
  public :: gfs_sigmaopen, gfs_sigmawrite, gfs_sigma2grid, gfs_grid2sigma

contains

  subroutine gfs_sigmaopen(unit,fname,head,data)
    !
    ! DESCRIPTION
    !  open GFS sigma file and read header & data
    !
    ! ARGUMENTS
    !  INPUT:
    !    unit  : file unit number
    !    fname : file name
    !  OUTPUT:
    !    head  : gfs sigio header (type:sigio_head)
    !    data  : gfs sigio data   (type:sigio_data)
    !
    !  LIBRARY
    !   sigio_r_module - sigio_rrohdc
    !
    integer(kind=i4b), intent(in)  :: unit
    character(len=*),  intent(in)  :: fname
    type(sigio_head),  intent(out) :: head
    type(sigio_data),  intent(out) :: data
    integer(kind=i4b) :: ierr

    call sigio_rrohdc(unit, fname, head, data, ierr)
    if (ierr /= 0) then
       print*, "file open error. ierr=", ierr
       call abort
    end if
    
  end subroutine gfs_sigmaopen

  subroutine gfs_sigmawrite(unit,fname,head,data)
    !
    ! DESCRIPTION
    !  open GFS sigma file and write header & data
    !
    ! ARGUMENTS
    !  INPUT:
    !    unit  : file unit number
    !    fname : file name
    !  OUTPUT:
    !    head  : gfs sigio header (type:sigio_head)
    !    data  : gfs sigio data   (type:sigio_data)
    !
    !  LIBRARY
    !   sigio_r_module - sigio_rwohdc
    !
    integer(kind=i4b), intent(in)     :: unit
    character(len=*),  intent(in)     :: fname
    type(sigio_head),  intent(inout)  :: head
    type(sigio_data),  intent(inout)  :: data
    integer(kind=i4b) :: ierr

    call sigio_rwohdc(unit, fname, head, data, ierr)
    if (ierr /= 0) then
       print*, "ERROR: SIGIO_RWOHDC ierr=", ierr
       call abort
    end if
    
  end subroutine gfs_sigmawrite

  subroutine gfs_sigma2grid(imax,jmax,kmax,maxwv,idrt,head,data,lon,lat,hs,ps,t,d,z,u,v,q,o3,cwc,yrev)
    !
    ! DESCRIPTION
    !  read spectral data of sigma file data and convert grid data
    !
    ! ARGUMENT
    !  INPUT:
    !    head  : gfs sigio header (type:sigio_head)
    !    data  : gfs sigio data   (type:sigio_data)
    !    idrt  : grid identifer
    !             idrt=0 - equally-spaced grid
    !             idrt=4 - gaussian grid
    !    maxwv : truncation max wave number
    !    yrev  : meridional grid direction flag (optional)
    !             True - N->S  False(default) - S->N
    !  OUTPUT:
    !    lon(imax)                : longitude
    !    lat(jmax)                : latitude
    !    hs(imax,jmax)            : geometric height [m]
    !    ps(imax,jmax)            : surface pressure [Pa]
    !    t(imax,jmax,kmax)        : dry temperature [K]
    !    d(imax,jmax,kmax)        : divergence
    !    z(imax,jmax,kmax)        : vorticity
    !    u(imax,jmax,kmax)        : horizotal wind [m/s]
    !    v(imax,jmax,kmax)        : meridional wind [m/s]
    !    q(imax,jmax,kmax)        : specific humidity
    !    o3(imax,jmax,kmax)       : ozone mixing ratio
    !    cwc(imax,jmax,kmax)      : cloud water mixing ratio
    !
    integer(kind=i4b), intent(in)  :: imax,jmax,kmax,maxwv,idrt
    type(sigio_head),  intent(in)  :: head
    type(sigio_data),  intent(in)  :: data
    real(kind=sp),     intent(out) :: lon(imax),lat(jmax)
    logical, intent(in), optional  :: yrev
    real(kind=sp),     intent(out) :: hs(imax,jmax),ps(imax,jmax)
    real(kind=sp),     intent(out) :: t(imax,jmax,kmax)
    real(kind=sp),     intent(out) :: d(imax,jmax,kmax),z(imax,jmax,kmax)
    real(kind=sp),     intent(out) :: u(imax,jmax,kmax),v(imax,jmax,kmax)
    real(kind=sp),     intent(out) :: q(imax,jmax,kmax),o3(imax,jmax,kmax),cwc(imax,jmax,kmax)

    integer(kind=i4b) :: i,k,j,ierr,ntrac
    real(kind=sp)     :: dd
    real(kind=sp) :: wlat(jmax)
    real(kind=sp), allocatable :: work(:,:), trac(:,:,:,:)
    logical :: isyrev

    ntrac = head%ntrac 
    allocate( trac(imax,jmax,kmax,ntrac) )
    if (present(yrev)) then
       isyrev = yrev
    else
       isyrev = .false.
    end if

    ! longitude
    dd = 360./real(imax)
    do i = 1, imax
       lon(i) = dd*(i-1)
    end do
    ! latitude
    if (idrt == 0) then
       d = 180./real(jmax)
       if (yrev) then
          do i = 1, jmax
             lat(i) = 90. - dd*(i-1)
          end do
       else
          do i = 1, jmax
             lat(i) = -90. + dd*(i-1)
          end do
       end if
    else if (idrt == 4) then
       call splat(4,jmax,lat,wlat)
       lat = asin(lat)*180./PI
       if (.not. isyrev) lat = lat(jmax:1:-1)
    else
       print*, "idrt must be 0 (equally-spaced grid including poles) or 4 (gaussian grid)"
       call abort
    end if
    
    call sptez(0,maxwv,idrt,imax,jmax,data%hs,hs,1)    ! geometric height
    call sptez(0,maxwv,idrt,imax,jmax,data%ps,ps,1)    ! surface pressure
    ps = exp(ps)*1.e3  !ln(KPa)->Pa  
    call sptezm(0,maxwv,idrt,imax,jmax,kmax,data%t,t,1) ! virtual temperature
    call sptezm(0,maxwv,idrt,imax,jmax,kmax,data%d,d,1) ! divergence
    call sptezm(0,maxwv,idrt,imax,jmax,kmax,data%z,z,1) ! vorticity  
    call sptezmv(0,maxwv,idrt,imax,jmax,kmax,data%d,data%z,u,v,1) ! zonal & meridional wind
    call sptezm(0,maxwv,idrt,imax,jmax,kmax*ntrac,data%q,trac,1) 
    !vertual temperature -> dry temperature
    call sigio_cnvtdv(imax*jmax,imax*jmax,kmax,head%idvc,head%idvm,head%ntrac, &
         &            ierr,t,trac,head%cpi,1)
    q   = trac(:,:,:,1)   ! specific humidity
    o3  = trac(:,:,:,2)   ! ozone
    cwc = trac(:,:,:,3)   ! cloud water

    if (.not. isyrev) then
       allocate( work(imax,jmax) )
       work = hs
       do j = 1, jmax
          hs(:,j) = work(:,jmax-j+1)
       end do
       work = ps
       do j = 1, jmax
          ps(:,j) = work(:,jmax-j+1)
       end do
       do k = 1, kmax
          work = t(:,:,k)
          do j = 1, jmax
             t(:,j,k) = work(:,jmax-j+1)
          end do
       end do
       do k = 1, kmax
          work = d(:,:,k)
          do j = 1, jmax
             d(:,j,k) = work(:,jmax-j+1)
          end do
       end do
       do k = 1, kmax
          work = z(:,:,k)
          do j = 1, jmax
             z(:,j,k) = work(:,jmax-j+1)
          end do
       end do
       do k = 1, kmax
          work = u(:,:,k)
          do j = 1, jmax
             u(:,j,k) = work(:,jmax-j+1)
          end do
       end do
       do k = 1, kmax
          work = v(:,:,k)
          do j = 1, jmax
             v(:,j,k) = work(:,jmax-j+1)
          end do
       end do
       do k = 1, kmax
          work = q(:,:,k)
          do j = 1, jmax
             q(:,j,k) = work(:,jmax-j+1)
          end do
       end do
       do k = 1, kmax
          work = o3(:,:,k)
          do j = 1, jmax
             o3(:,j,k) = work(:,jmax-j+1)
          end do
       end do
       do k = 1, kmax
          work = cwc(:,:,k)
          do j = 1, jmax
             cwc(:,j,k) = work(:,jmax-j+1)
          end do
       end do
    end if

  end subroutine gfs_sigma2grid

  subroutine gfs_grid2sigma(imax,jmax,kmax,maxwv,idrt,head,data,hs,ps,t,d,z,q,o3,cwc,yrev)
    !
    ! DESCRIPTION
    !  read spectral data of sigma file data and convert grid data
    !
    ! ARGUMENT
    !  INPUT:
    !    head  : gfs sigio header (type:sigio_head)
    !    data  : gfs sigio data   (type:sigio_data)
    !    idrt  : grid identifer
    !             idrt=0 - equally-spaced grid
    !             idrt=4 - gaussian grid
    !    maxwv : truncation max wave number
    !    yrev  : meridional grid direction flag (optional)
    !             True - N->S  False(default) - S->N
    !  INPUT:
    !    hs(imax,jmax)       : geometric height [m]
    !    ps(imax,jmax)       : surface pressure [Pa]
    !    t(imax,jmax,kmax)   : dry temperature [K]
    !    d(imax,jmax,kmax)   : divergence
    !    z(imax,jmax,kmax)   : vorticity
    !    q(imax,jmax,kmax)   : specific humidity
    !    o3(imax,jmax,kmax)  : ozone mixing ratio
    !    cwc(imax,jmax,kmax) : cloud water mixing ratio
    !
    integer(kind=i4b), intent(in)  :: imax,jmax,kmax,maxwv,idrt
    type(sigio_head),  intent(in)  :: head
    type(sigio_data),  intent(out) :: data
    logical, intent(in), optional  :: yrev
    real(kind=sp),     intent(out) :: hs(imax,jmax),ps(imax,jmax)
    real(kind=sp),     intent(out) :: t(imax,jmax,kmax)
    real(kind=sp),     intent(out) :: d(imax,jmax,kmax),z(imax,jmax,kmax)
    real(kind=sp),     intent(out) :: q(imax,jmax,kmax),o3(imax,jmax,kmax),cwc(imax,jmax,kmax)

    integer(kind=i4b) :: i,k,j,ierr,ntrac
    real(kind=sp)     :: dd
    real(kind=sp), allocatable :: work(:,:), trac(:,:,:,:)
    logical :: isyrev

    ! deallocate data
    call sigio_axdata(data,ierr)
    ! allocate data
    call sigio_aldata(head,data,ierr)
    if ( ierr/= 0 ) then
       print*, "err, sigio_aldata. ierr=",ierr
       call abort
    end if

    ntrac = head%ntrac 
    allocate( trac(imax,jmax,kmax,ntrac) )
    if (present(yrev)) then
       isyrev = yrev
    else
       isyrev = .false.
    end if

    if (.not. isyrev) then
       allocate( work(imax,jmax) )
       work = hs
       do j = 1, jmax
          hs(:,j) = work(:,jmax-j+1)
       end do
       work = ps
       do j = 1, jmax
          ps(:,j) = work(:,jmax-j+1)
       end do
       do k = 1, kmax
          work = t(:,:,k)
          do j = 1, jmax
             t(:,j,k) = work(:,jmax-j+1)
          end do
       end do
       do k = 1, kmax
          work = d(:,:,k)
          do j = 1, jmax
             d(:,j,k) = work(:,jmax-j+1)
          end do
       end do
       do k = 1, kmax
          work = z(:,:,k)
          do j = 1, jmax
             z(:,j,k) = work(:,jmax-j+1)
          end do
       end do
       do k = 1, kmax
          work = q(:,:,k)
          do j = 1, jmax
             q(:,j,k) = work(:,jmax-j+1)
          end do
       end do
       do k = 1, kmax
          work = o3(:,:,k)
          do j = 1, jmax
             o3(:,j,k) = work(:,jmax-j+1)
          end do
       end do
       do k = 1, kmax
          work = cwc(:,:,k)
          do j = 1, jmax
             cwc(:,j,k) = work(:,jmax-j+1)
          end do
       end do
    end if
    
    call sptez(0,maxwv,idrt,imax,jmax,data%hs,hs,-1)    ! geometric height
    ps = log(ps*1.e-3)     !Pa->ln(KPa)
    call sptez(0,maxwv,idrt,imax,jmax,data%ps,ps,-1)    ! surface pressure
    trac(:,:,:,1) = q                                   ! specific humidity
    trac(:,:,:,2) = o3                                  ! ozone
    trac(:,:,:,3) = cwc                                 ! cloud water
    !dry temperature -> virtual temperature
    call sigio_cnvtdv(imax*jmax,imax*jmax,kmax,head%idvc,head%idvm,head%ntrac, &
         &            ierr,t,trac,head%cpi,-1)
    call sptezm(0,maxwv,idrt,imax,jmax,kmax,data%t,t,-1) ! virtual temperature
    call sptezm(0,maxwv,idrt,imax,jmax,kmax,data%d,d,-1) ! divergence
    call sptezm(0,maxwv,idrt,imax,jmax,kmax,data%z,z,-1) ! vorticity  
    call sptezm(0,maxwv,idrt,imax,jmax,kmax*ntrac,data%q,trac,-1) 

  end subroutine gfs_grid2sigma


end module gfs_module
