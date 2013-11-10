program main
  use type_module
  use constant_module, only : PI=>math_pi  
  use sigio_module
  use gfs_module
  use p2p_module
  use calcmet_module
  use ip_module
  use vortex_module
  !for debug
  use debug_module
  implicit none
  !!
  ! PROGRAM vortex_move
  ! 
  ! DESCRIPTION
  !  This program read gfs sigma file spectral data and move 
  ! the hurricane vortex component to spcifield position
  ! and write gfs sigma file.
  !
  ! FILES
  !  INPUT FILE:
  !   unit 11: namelist file
  !   unit 12: input gfs sigma file
  !   unit 21: output gfs sigma file
  !
  ! LIBRARY
  !  sigio - sigio_modpr
  !
  !!
  ! PARAMETER - irmax,jrmax,dxyr
  !  size of subdomain grid to applay vortex removal scheme
  !  41x41 horizontal grids with grid intervals of 1 degree (Liu et al. (2000))
  !
  ! AUTHOR
  !  Tetsuro Miyachi
  !   - https://bitbucket.org/tmiyachi
  !   - https://github.com/tmiyachi
  !   
  !!
  integer(kind=i4b), parameter :: irmax=41,jrmax=41,ijrmax=irmax*jrmax
  real(kind=sp),     parameter :: dxyr=1.

  integer(kind=i4b), parameter :: fnm=11,fni=12,fno=21

  integer(kind=i4b) :: imax,jmax,kmax,nmax,irgmax,jrgmax
  integer(kind=i4b) :: maxwv,ntrac,nvcoord,idsl,idvc
  real(kind=sp),allocatable :: vcoord(:,:) 
  type(sigio_head) :: headi,heado
  type(sigio_data) :: datai,datao

  real(kind=sp) :: clon_obs,clat_obs,clon_new,clat_new,clon,clat,dlon,dlat

  real(kind=sp), allocatable :: lon(:),lat(:)
  real(kind=sp), allocatable :: hsg(:,:),psg(:,:),tg(:,:,:)
  real(kind=sp), allocatable :: dg(:,:,:),zg(:,:,:),ug(:,:,:),vg(:,:,:)
  real(kind=sp), allocatable :: qg(:,:,:),o3g(:,:,:),cwcg(:,:,:)

  real(kind=sp)              :: lonr(irmax),latr(jrmax)
  real(kind=sp), allocatable :: hsr(:,:),datrm(:,:,:),datrp(:,:,:)
  real(kind=sp), allocatable :: prm(:,:,:),psref(:,:),prp(:,:,:)

  real(kind=sp), allocatable :: lu(:,:),lv(:,:),env(:,:,:),vrtex(:,:,:)

  integer(kind=i4b) :: klu,klv,ix1,ix2,jy1,jy2,ixc,jyc

  integer(kind=i4b) :: i,j,k,n,ierr,dk
  integer(kind=i4b) :: nv,nvmax,nvps,nvt,nvu,nvv,nvq

  namelist/tcinfo/clon_obs,clat_obs
  namelist/param/clon_new,clat_new

  !!
  !! 0) read namelist
  !!
  clon_new = -9999
  clat_new = -9999
  write(*,*) "[0] READ NAMELIST"
  open(fnm,file="vortex_move_namelist",iostat=ierr)
  if ( ierr /= 0 ) then
     print*, "namelist file open error. ierr=", ierr
     call abort
  end if
  read(fnm,nml=tcinfo)
  read(fnm,nml=param)
  close(fnm)
  write(*,'("  observed TC center clon,clat=  ",f5.1,f5.1)') clon_obs,clat_obs
  write(*,*)

  !!
  !! 1) read global grid data
  !!
  write(*,*) "[1] OPEN INPUT GFS SIGMA FILE"

  call gfs_sigmaopen(fno,'input.sigma',headi,datai)
  imax    = headi%lonb
  jmax    = headi%latb
  kmax    = headi%levs
  ntrac   = headi%ntrac
  maxwv   = headi%jcap
  nvcoord = headi%nvcoord
  idvc    = headi%idvc
  idsl    = headi%idsl
  allocate( vcoord(kmax+1,nvcoord) )
  vcoord  = headi%vcoord
  
  allocate( lon(imax), lat(jmax) )
  allocate( hsg(imax,jmax), psg(imax,jmax) )
  allocate( tg(imax,jmax,kmax) )
  allocate( dg(imax,jmax,kmax), zg(imax,jmax,kmax) )
  allocate( ug(imax,jmax,kmax), vg(imax,jmax,kmax) )
  allocate( qg(imax,jmax,kmax), o3g(imax,jmax,kmax), cwcg(imax,jmax,kmax) )

  call gfs_sigma2grid(imax,jmax,kmax,maxwv,4,headi,datai,lon,lat,hsg,psg,tg, &
       &            dg,zg,ug,vg,qg,o3g,cwcg,yrev=.false.)
  write(*,'("  imax,jmax,kmax= ",i5,i5,i5)') imax,jmax,kmax
  write(*,'("  maxwv= ",i0)') maxwv
  write(*,*) 


  !!
  !! 2) interpolate to the regional subdomain 
  !!
  write(*,*) "[2] INTERPOLATION TO THE REGIONAL SUBDOMAIN"

  ixc = searchidx(lon,clon_obs,0)    
  do i = 1, irmax
     lonr(i) = lon(ixc) + dxyr*real(-irmax/2 + (i-1))
  end do
  jyc = searchidx(lat,clat_obs,0)    
  do j = 1, jrmax
     latr(j) = lat(jyc) + dxyr*real(-jrmax/2 + (j-1))
  end do

  nvmax   = 1 + 5*kmax + 3*kmax
  dk = kmax - 1
  allocate( hsr(irmax,jrmax))
  allocate( datrm(irmax,jrmax,nvmax), datrp(irmax,jrmax,nvmax))
  call interp2d(imax,jmax,1,lon,lat,hsg,irmax,jrmax,lonr,latr,hsr)
  nv   = 1
  nvps = nv
  call interp2d(imax,jmax,1,lon,lat,psg,irmax,jrmax,lonr,latr,datrm(:,:,nv))
  nv  = nv + 1
  nvt = nv
  call interp2d(imax,jmax,kmax,lon,lat,tg,irmax,jrmax,lonr,latr,datrm(:,:,nv:nv+dk))
  nv = nv + kmax
  call interp2d(imax,jmax,kmax,lon,lat,dg,irmax,jrmax,lonr,latr,datrm(:,:,nv:nv+dk))
  nv = nv + kmax
  call interp2d(imax,jmax,kmax,lon,lat,zg,irmax,jrmax,lonr,latr,datrm(:,:,nv:nv+dk))
  nv  = nv + kmax
  nvu = nv 
  call interp2d(imax,jmax,kmax,lon,lat,ug,irmax,jrmax,lonr,latr,datrm(:,:,nv:nv+dk))
  nv  = nv + kmax
  nvv = nv 
  call interp2d(imax,jmax,kmax,lon,lat,vg,irmax,jrmax,lonr,latr,datrm(:,:,nv:nv+dk))
  nv  = nv + kmax
  nvq = nv
  call interp2d(imax,jmax,kmax,lon,lat,qg,irmax,jrmax,lonr,latr,datrm(:,:,nv:nv+dk))
  nv = nv + kmax
  call interp2d(imax,jmax,kmax,lon,lat,o3g,irmax,jrmax,lonr,latr,datrm(:,:,nv:nv+dk))
  nv = nv + kmax
  call interp2d(imax,jmax,kmax,lon,lat,cwcg,irmax,jrmax,lonr,latr,datrm(:,:,nv:nv+dk))

  write(*,'(i4" x ",i0," domain ",f4.1," degrees interval")'),irmax,jrmax,dxyr
  write(*,'("  domain center is     ",f7.2,f7.2)'),lon(ixc),lat(jyc)
  write(*,'("  lonr(1),lonr(irmax)= ",f7.2,f7.2)'),lonr(1),lonr(irmax)
  write(*,'("  latr(1),latr(jrmax)= ",f7.2,f7.2)'),latr(1),latr(jrmax)
  write(*,*)
  

  !!
  !! 3) vertical interpolate to the constant pressure level
  !!
  write(*,*) "[3] VERTICAL INTERPOLATION TO THE CONSTANT PRESSURE LEVEL"

  !! 3.1) calculate input model pressure level
  allocate( prm(irmax,jrmax,kmax), psref(irmax,jrmax), prp(irmax,jrmax,kmax) )
  call sigio_modpr(ijrmax,ijrmax,kmax,nvcoord,idvc,idsl,vcoord,ierr,ps=datrm(:,:,nvps),pm=prm)

  !! 3.2) calculate constant pressure level for vertical interpolation
  psref = maxval(datrm(:,:,nvps))    ! referrence pressure
  call sigio_modpr(ijrmax,ijrmax,kmax,nvcoord,idvc,idsl,vcoord,ierr,ps=psref,pm=prp)
  
  ! for debug  
  call write_grads(21,'check_mlev_total.bin',irmax,jrmax,kmax,1,8,lonr,latr,datrm)

  !! 3.3) vertically interpolate data from model levels to constant pressure levels
  datrp(:,:,nvps) = datrm(:,:,nvps)
  do nv = nvt, nvmax, kmax
     if (nv == nvt) then
        call p2p_extrapolate_T(ijrmax,kmax,datrm(:,:,nv:nv+dk),prm,datrm(:,:,nvps),   &
             &                 kmax,prp,datrp(:,:,nvps),hsr,datrp(:,:,nv:nv+dk))
     else if (nv == nvq) then
        datrm(:,:,nv:nv+dk) = calcmet_rh(datrm(:,:,nvt:nvt+dk),datrm(:,:,nv:nv+dk),prm)   !SPH->RH        
        call p2p(ijrmax,kmax,datrm(:,:,nv:nv+dk),prm,datrm(:,:,nvps),kmax,prp,datrp(:,:,nv:nv+dk))
        datrp(:,:,nv:nv+dk) = min(datrp(:,:,nv:nv+dk),100.)
        datrp(:,:,nv:nv+dk) = max(datrp(:,:,nv:nv+dk),0.)
        datrp(:,:,nv:nv+dk) = calcmet_q(datrp(:,:,nvt:nvt+dk),datrp(:,:,nv:nv+dk),prp)    !RH->SPH
     else        
        call p2p(ijrmax,kmax,datrm(:,:,nv:nv+dk),prm,datrm(:,:,nvps),kmax,prp,datrp(:,:,nv:nv+dk))
     end if
  end do
  ! convert Ps->Msl
  !datrp(:,:,nvps) = calcmet_msl(datrm(:,:,nvt:nvt+dk),hsr,datrp(:,:,nvps),prm)
  write(*,'("  refference pressure is ",f9.2)') psref(1,1) 
  write(*,'("  prm(1,1,1)= ",f9.2,"  prp(1,1,1)= ",f9.2)') prm(1,1,1),prp(1,1,1)
  write(*,*)
  ! for debug
  !  call write_grads(21,'check_plev_total.bin',irmax,jrmax,kmax,1,8,lonr,latr,datrp)


  !!
  !! 4) separate vortex field
  !!
  write(*,*) "[4] COMPUTE VORTEX FIELD AND ENVIRONMENTAL FIELD"
  k = searchidx(prp(1,1,:),85000.,0)  
  allocate( lu(irmax,jrmax), lv(irmax,jrmax) )
  allocate( env(irmax,jrmax,nvmax),vrtex(irmax,jrmax,nvmax))
  lu = datrp(:,:,nvu+k-1)
  lv = datrp(:,:,nvv+k-1)
  write(*,'("  determine filter domain at ",f9.2,"Pa, k=",i0)'),prp(1,1,k),k 
  write(*,*)
  call separate_env_vortex(irmax,jrmax,nvmax,lonr,latr,datrp,lu,lv,clon_obs,clat_obs,env,vrtex,clon,clat)
  write(*,*)
  ! for debug
  !  call write_grads(21,'check_plev_total.bin',irmax,jrmax,kmax,1,8,lonr,latr,datrp)
  !  call write_grads(21,'check_plev_env.bin',irmax,jrmax,kmax,1,8,lonr,latr,env)
  !  call write_grads(21,'check_plev_vrtex.bin',irmax,jrmax,kmax,1,8,lonr,latr,vrtex)

  !!
  !! 5) vortex relocation
  !!
  write(*,*) "[5] VORTEX RELOCATION"
  if (clon_new==-9999) clon_new = clon_obs
  if (clat_new==-9999) clat_new = clat_obs
  
  dlon = clon_new - clon
  dlat = clat_new - clat
  call interp2d(irmax,jrmax,nvmax,lonr,latr,vrtex,irmax,jrmax,lonr-dlon,latr-dlat,datrp,undef_value=0.)
  datrp = env + datrp 
  write(*,'("  vortex TC center clon,clat=           ",f6.1,f6.1)') clon,clat
  write(*,'("  relocated TC center clon_new,clat_new=",f6.1,f6.1)') clon_new,clat_new
  write(*,'("  difference of TC center dlon,dlat=    ",f6.2,f6.2)') dlon,dlat
  write(*,*)
  ! for debug
  !  call write_grads(21,'check_plev_relocated.bin',irmax,jrmax,kmax,1,8,lonr,latr,datrp)


  !!
  !! 6) vertical interpolation to model level
  !!
  write(*,*) "[6] VERTICAL INTERPOLATION TO THE MODEL LEVEL"  
  write(*,*)
  datrm(:,:,nvps) = datrp(:,:,nvps)
  do nv = nvt, nvmax, kmax
     if (nv == nvt) then
        call p2p_extrapolate_T(ijrmax,kmax,datrp(:,:,nv:nv+dk),prp,datrp(:,:,nvps),    &
             &                 kmax,prm,datrm(:,:,nvps),hsr,datrm(:,:,nv:nv+dk))
     else if (nv == nvq) then
        datrp(:,:,nv:nv+dk) = calcmet_rh(datrp(:,:,nvt:nvt+dk),datrp(:,:,nv:nv+dk),prp)   !SPH->RH        
        call p2p(ijrmax,kmax,datrp(:,:,nv:nv+dk),prp,datrp(:,:,nvps),kmax,prm,datrm(:,:,nv:nv+dk))
        datrm(:,:,nv:nv+dk) = calcmet_q(datrm(:,:,nvt:nvt+dk),datrm(:,:,nv:nv+dk),prm)    !RH->SPH
        datrm(:,:,nv:nv+dk) = min(datrm(:,:,nv:nv+dk),100.)
        datrm(:,:,nv:nv+dk) = max(datrm(:,:,nv:nv+dk),0.)
     else        
        call p2p(ijrmax,kmax,datrp(:,:,nv:nv+dk),prp,datrp(:,:,nvps),kmax,prm,datrm(:,:,nv:nv+dk))        
     end if
  end do
  ! convert MSL->Ps
  !datrm(:,:,nvps) = calcmet_ps(datrm(:,:,nvt:nvt+dk),hsr,datrm(:,:,nvps),prm)
  ! for debug
  !call write_grads(21,'check_mlev_relocated.bin',irmax,jrmax,kmax,1,8,lonr,latr,datrm)


  !!
  !! 7) interpolation to global grid
  !!
  write(*,*) "[7] INTERPOLATION TO GLOBAL GRID"
  !                                     ! R(1,jrmax)        -            R(irmax,jrmax)
  ix1 = searchidx(lon,lonr(1),-1)       !        G(ix1,jy2) - G(ix2,iy2)
  ix2 = searchidx(lon,lonr(irmax),1)    !    |        |           |         |
  jy1 = searchidx(lat,latr(1),-1)       !        G(ix1,jy1) - G(ix2,jy1)
  jy2 = searchidx(lat,latr(jrmax),1)    ! R(1,1)            -            R(irmax,1)
  irgmax = ix2-ix1+1
  jrgmax = jy2-jy1+1

  write(*,'("  lon(ix1),lon(ix2),lat(jy1),lag(jy2)= ",f7.2,f7.2,f6.2,f6.2)'),lon(ix1),lon(ix2),lat(jy1),lat(jy2)
  write(*,*) 
  nv = 1  
  call interp2d(irmax,jrmax,1,lonr,latr,datrm(:,:,nv),            &
       &        irgmax,jrgmax,lon(ix1:ix2),lat(jy1:jy2),psg(ix1:ix2,jy1:jy2))  
  nv = nv + 1
  call interp2d(irmax,jrmax,kmax,lonr,latr,datrm(:,:,nv:nv+dk),   &
       &        irgmax,jrgmax,lon(ix1:ix2),lat(jy1:jy2),tg(ix1:ix2,jy1:jy2,:))
  nv = nv + kmax
  call interp2d(irmax,jrmax,kmax,lonr,latr,datrm(:,:,nv:nv+dk),   &
       &        irgmax,jrgmax,lon(ix1:ix2),lat(jy1:jy2),dg(ix1:ix2,jy1:jy2,:))
  nv = nv + kmax
  call interp2d(irmax,jrmax,kmax,lonr,latr,datrm(:,:,nv:nv+dk),   &
       &        irgmax,jrgmax,lon(ix1:ix2),lat(jy1:jy2),zg(ix1:ix2,jy1:jy2,:))
  nv = nv + kmax
  call interp2d(irmax,jrmax,kmax,lonr,latr,datrm(:,:,nv:nv+dk),   &
       &        irgmax,jrgmax,lon(ix1:ix2),lat(jy1:jy2),ug(ix1:ix2,jy1:jy2,:))
  nv = nv + kmax
  call interp2d(irmax,jrmax,kmax,lonr,latr,datrm(:,:,nv:nv+dk),   &
       &        irgmax,jrgmax,lon(ix1:ix2),lat(jy1:jy2),vg(ix1:ix2,jy1:jy2,:))
  nv = nv + kmax
  call interp2d(irmax,jrmax,kmax,lonr,latr,datrm(:,:,nv:nv+dk),   &
       &        irgmax,jrgmax,lon(ix1:ix2),lat(jy1:jy2),qg(ix1:ix2,jy1:jy2,:))
  nv = nv + kmax
  call interp2d(irmax,jrmax,kmax,lonr,latr,datrm(:,:,nv:nv+dk),   &
       &        irgmax,jrgmax,lon(ix1:ix2),lat(jy1:jy2),o3g(ix1:ix2,jy1:jy2,:))
  nv = nv + kmax
  call interp2d(irmax,jrmax,kmax,lonr,latr,datrm(:,:,nv:nv+dk),   &
       &        irgmax,jrgmax,lon(ix1:ix2),lat(jy1:jy2),cwcg(ix1:ix2,jy1:jy2,:))


  !!
  !! 8) write output gfs sigma file
  !!
  write(*,*) "WRITE OUTPUT GFS SIGMA FILE"
  call gfs_grid2sigma(imax,jmax,kmax,maxwv,4,headi,datao,hsg,psg,tg,dg,zg,qg,o3g,cwcg,yrev=.false.)
  call gfs_sigmawrite(fno,'output.sigma',headi,datao)


end program main




