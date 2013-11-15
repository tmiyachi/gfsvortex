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
  !use debug_module
  implicit none
  !!
  ! PROGRAM vortex_move
  ! 
  ! DESCRIPTION
  !  This program split the field on the GFS model level into the hurricane vortex component
  !  and the environtal component, and move the vortex componet to the new location.
  !  This program contains following steps.
  !
  !   1. convert surface pressure to sea level pressure to avoid the steep montain effects.
  !   2. vertically interpolate from model level to constant pressure level
  !      within regional subdomain.
  !   3. horizontally interpolate the regional field from the gaussian grid to the equally
  !      spaced lat-lon grid on the each constant pressure level.
  !   4. split the total field into the environmental field and the vortex field on the each constant
  !      pressure level using the scheme proposed by Kurihara et al.(1992,1995)
  !   5. relocate the vortex field to the new TC center.
  !   6. horizontally interpolate the the regional field from the equally spaced lat-lon grid
  !      to the gaussian grid on the each constant pressure level.
  !   7. convert sea level pressure to surface pressure.
  !   8. vertically interpolate from the constant pressure level to the model level coressponding 
  !      to the relocated surface pressure.
  !
  ! FILES
  !  INPUT FILE:
  !   unit 11: namelist file
  !   unit 12: input gfs sigma file
  !   unit 21: output gfs sigma file
  !
  ! NAMELIST
  !  tcinfo
  !    clon_obs,clat_obs : observed TC center to use first-guess TC center position 
  !  param
  !    clon_new,clat_new : new relocated TC center position. observed TC center is used in default
  !
  ! LIBRARY
  !  sigio - sigio_modpr
  !
  ! AUTHOR
  !  Tetsuro Miyachi
  !   - https://bitbucket.org/tmiyachi
  !   - https://github.com/tmiyachi
  !
  ! HISTORY
  !  2013.11.07 - initial version
  !  see the commit log of the repositories
  !   
  !!
  ! PARAMETER - irmax,jrmax,dxyr
  !  size of subdomain grid to applay vortex removal scheme
  !  41x41 horizontal grids with grid intervals of 1 degree (Liu et al. (2000))
  !!
  integer(kind=i4b), parameter :: fnm=11,fni=12,fno=21

  integer(kind=i4b), parameter :: irmax=41,jrmax=41
  real(kind=sp),     parameter :: dxyr=1.

  integer(kind=i4b)          :: imax,jmax,kmax
  integer(kind=i4b)          :: maxwv,ntrac,nvcoord,idsl,idvc,idate(4)
  real(kind=sp),allocatable  :: vcoord(:,:) 
  type(sigio_head)           :: headi,heado
  type(sigio_data)           :: datai,datao

  real(kind=sp) :: clon_obs,clat_obs,clon_new,clat_new

  real(kind=sp), allocatable :: lon(:),lat(:)
  real(kind=sp), allocatable :: hsg(:,:),psg(:,:),tg(:,:,:)
  real(kind=sp), allocatable :: dg(:,:,:),zg(:,:,:),ug(:,:,:),vg(:,:,:)
  real(kind=sp), allocatable :: qg(:,:,:),o3g(:,:,:),cwcg(:,:,:)

  real(kind=sp)              :: psref(1)
  real(kind=sp), allocatable :: plev(:)

  integer(kind=i4b)          :: irgmax,jrgmax,ijrgmax
  real(kind=sp), allocatable :: datrgp(:,:,:), prgm(:,:,:), prgp(:,:,:)

  real(kind=sp)              :: lonr(irmax),latr(jrmax)
  real(kind=sp), allocatable :: datrp(:,:,:)
  real(kind=sp), allocatable :: lu(:,:),lv(:,:),env(:,:,:),vrtex(:,:,:)
  
  real(kind=sp) :: clon,clat,dlon,dlat
  

  integer(kind=i4b) :: klu,klv,ix1,ix2,jy1,jy2,ixc,jyc

  integer(kind=i4b) :: i,j,k,n,ierr,dk
  integer(kind=i4b) :: nv,nvmax,nvt,nvu,nvv,nvq,nvslp

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
  write(*,'("  new TC center clon,clat=       ",f5.1,f5.1)') clon_new,clat_new
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
  idate   = headi%idate
  allocate( vcoord(kmax+1,nvcoord) )
  vcoord  = headi%vcoord
  nvmax   = 1 + 5*kmax + 3*kmax         !ps(slp),T,d,z,u,v,q,o3,cwc
  dk = kmax - 1
  
  allocate( lon(imax), lat(jmax) )
  allocate( hsg(imax,jmax), psg(imax,jmax) )
  allocate( tg(imax,jmax,kmax) )
  allocate( dg(imax,jmax,kmax), zg(imax,jmax,kmax) )
  allocate( ug(imax,jmax,kmax), vg(imax,jmax,kmax) )
  allocate( qg(imax,jmax,kmax), o3g(imax,jmax,kmax), cwcg(imax,jmax,kmax) )

  call gfs_sigma2grid(imax,jmax,kmax,maxwv,4,headi,datai,lon,lat,hsg,psg,tg, &
       &            dg,zg,ug,vg,qg,o3g,cwcg,yrev=.false.)
  write(*,'("  imax,jmax,kmax= ",i5,i5,i5)') imax,jmax,kmax
  write(*,'("  maxwv= ",i0, "  idate=",i5,i3,i3,i3)') maxwv,idate(4),idate(2),idate(3),idate(1)
  write(*,*) 


  !!
  !! 2) calculate regional subdomain
  !!
  write(*,*) "[2] CALCULATE REGIONAL SUBDOMAIN"
  ixc = searchidx(lon,clon_obs,0)    
  do i = 1, irmax
     lonr(i) = lon(ixc) + dxyr*real(-irmax/2 + (i-1))
  end do
  jyc = searchidx(lat,clat_obs,0)    
  do j = 1, jrmax
     latr(j) = lat(jyc) + dxyr*real(-jrmax/2 + (j-1))
  end do
  write(*,'(i4" x ",i0," domain ",f4.1," degrees interval")'),irmax,jrmax,dxyr
  write(*,'("  domain center is     ",f7.2,f7.2)'),lon(ixc),lat(jyc)
  write(*,'("  lonr(1),lonr(irmax)= ",f7.2,f7.2)'),lonr(1),lonr(irmax)
  write(*,'("  latr(1),latr(jrmax)= ",f7.2,f7.2)'),latr(1),latr(jrmax)
  write(*,*) 
  

  !!
  !! 3) vertically interpolate from the model level to the constant pressure level
  !!
  write(*,*) "[3] VERTICAL INTERPOLATION TO THE CONSTANT PRESSURE LEVEL"

  !!
  !! need two rows & columns outside reigonal subdomain to use quasi-cubic interpolation
  !!
  !                                       ! G(ix+1,jy2-1)    -         G(ix2-1,jy2-1)
  ix1 = searchidx(lon,lonr(1),1) - 1      !       R(1,jrmax) - R(irmax,jrmax)
  ix2 = searchidx(lon,lonr(irmax),-1) + 1 ! |     |          |                      |
  jy1 = searchidx(lat,latr(1),1) - 1      !       R(1,1)     - R(irmax,1)
  jy2 = searchidx(lat,latr(jrmax),-1) + 1 ! G(ix+1,jy1+1)    -         G(ix2-1,jy1+1)
  irgmax = ix2-ix1+1                         
  jrgmax = jy2-jy1+1                      
  ijrgmax = irgmax*jrgmax

  allocate( datrgp(irgmax,jrgmax,nvmax) )
  allocate( plev(kmax) )
  allocate( prgm(irgmax,jrgmax,kmax), prgp(irgmax,jrgmax,kmax) )

  !! 3.1) calculate input model pressure level
  call sigio_modpr(ijrgmax,ijrgmax,kmax,nvcoord,idvc,idsl,vcoord,ierr,ps=psg(ix1:ix2,jy1:jy2),pm=prgm)

  !! 3.2) calculate constant pressure level for vertical interpolation
  psref = maxval(psg(ix1:ix2,jy1:jy2))    ! referrence pressure
  call sigio_modpr(1,1,kmax,nvcoord,idvc,idsl,vcoord,ierr,ps=psref,pm=plev)
  do k = 1, kmax
     prgp(:,:,k) = plev(k)
  end do

  !! 3.3) convert surface pressure to sea level pressure
  nvslp = 1
  datrgp(:,:,nvslp) = calcmet_ps2slp(tg(ix1:ix2,jy1:jy2,:),hsg(ix1:ix2,jy1:jy2),psg(ix1:ix2,jy1:jy2),prgm) 

  !! 3.4) vertically interpolate data from model levels to constant pressure levels
  nv  = 2
  nvt = nv
  call p2p_extrapolate_T(ijrgmax,kmax,tg(ix1:ix2,jy1:jy2,:),prgm,psg(ix1:ix2,jy1:jy2),   &
       &                 kmax,prgp,psg(ix1:ix2,jy1:jy2),hsg(ix1:ix2,jy1:jy2),datrgp(:,:,nv:nv+dk))
  nv = nv + kmax
  call p2p(ijrgmax,kmax,dg(ix1:ix2,jy1:jy2,:),prgm,psg(ix1:ix2,jy1:jy2),kmax,prgp,datrgp(:,:,nv:nv+dk))
  nv = nv + kmax
  call p2p(ijrgmax,kmax,zg(ix1:ix2,jy1:jy2,:),prgm,psg(ix1:ix2,jy1:jy2),kmax,prgp,datrgp(:,:,nv:nv+dk))
  nv = nv + kmax
  nvu = nv
  call p2p(ijrgmax,kmax,ug(ix1:ix2,jy1:jy2,:),prgm,psg(ix1:ix2,jy1:jy2),kmax,prgp,datrgp(:,:,nv:nv+dk))
  nv = nv + kmax
  nvv = nv
  call p2p(ijrgmax,kmax,vg(ix1:ix2,jy1:jy2,:),prgm,psg(ix1:ix2,jy1:jy2),kmax,prgp,datrgp(:,:,nv:nv+dk))
  nv = nv + kmax
  nvq = nv  
  call p2p_q2rh(ijrgmax,kmax,qg(ix1:ix2,jy1:jy2,:),tg(ix1:ix2,jy1:jy2,:),prgm,psg(ix1:ix2,jy1:jy2),  &
       &        kmax,prgp,psg(ix1:ix2,jy1:jy2),hsg(ix1:ix2,jy1:jy2),datrgp(:,:,nv:nv+dk))  
  nv = nv + kmax
  call p2p(ijrgmax,kmax,o3g(ix1:ix2,jy1:jy2,:),prgm,psg(ix1:ix2,jy1:jy2),kmax,prgp,datrgp(:,:,nv:nv+dk))
  nv = nv + kmax
  call p2p(ijrgmax,kmax,cwcg(ix1:ix2,jy1:jy2,:),prgm,psg(ix1:ix2,jy1:jy2),kmax,prgp,datrgp(:,:,nv:nv+dk))

  write(*,'("  refference pressure is ",f9.2)') psref(1) 
  write(*,'("  prgm(1,1,1)=       ",f9.2,"  prgp(1,1,1)=     ",f9.2)') prgm(1,1,1),prgp(1,1,1)
  write(*,'("  datrgp(1,1,nvslp)= ",f9.2,"  datrgp(1,1,nvt)= ",f9.2)') datrgp(1,1,nvslp),datrgp(1,1,nvt)
  write(*,'("  datrgp(1,1,nvu)=   ",f9.2,"  datrgp(1,1,nvv)= ",f9.2)') datrgp(1,1,nvu),datrgp(1,1,nvv)
  write(*,'("  datrgp(1,1,nvq)=   ",f9.2)') datrgp(1,1,nvq)
  write(*,*)


  !!
  !! 4) horizontal interpolation
  !!
  write(*,*) "[4] HORIZONTAL INTERPOLATION FROM GAUSSIAN GRID TO REGIONAL GRID"
  write(*,*) 
  allocate( datrp(irmax,jrmax,nvmax) )
  call interp2d_quasicubic(irgmax,jrgmax,nvmax,lon(ix1:ix2),lat(jy1:jy2),datrgp,irmax,jrmax,lonr,latr,datrp)


  !!
  !! 5) separate vortex field
  !!
  write(*,*) "[5] COMPUTE VORTEX FIELD AND ENVIRONMENTAL FIELD"
  k = searchidx(prgp(1,1,:),85000.,0)  
  allocate( lu(irmax,jrmax), lv(irmax,jrmax) )
  allocate( env(irmax,jrmax,nvmax),vrtex(irmax,jrmax,nvmax))
  lu = datrp(:,:,nvu+k-1)
  lv = datrp(:,:,nvv+k-1)
  write(*,'("  determine filter domain at ",f9.2,"Pa, k=",i0)'),prgp(1,1,k),k 
  write(*,*)

  call separate_env_vortex(irmax,jrmax,nvmax,lonr,latr,datrp,lu,lv,clon_obs,clat_obs,env,vrtex,clon,clat)
  write(*,*)

  ! for debug
  !call write_grads(21,'check_plev_total.bin',irmax,jrmax,kmax,1,8,lonr,latr,datrp)
  !call write_grads(21,'check_plev_env.bin',irmax,jrmax,kmax,1,8,lonr,latr,env)
  !call write_grads(21,'check_plev_vrtex.bin',irmax,jrmax,kmax,1,8,lonr,latr,vrtex)


  !!
  !! 6) vortex relocation
  !!
  write(*,*) "[6] VORTEX RELOCATION"
  if (clon_new==-9999) clon_new = clon_obs
  if (clat_new==-9999) clat_new = clat_obs
  
  dlon = clon_new - clon
  dlat = clat_new - clat
  call interp2d_quasicubic(irmax,jrmax,nvmax,lonr,latr,vrtex,irmax,jrmax,& 
       &                   lonr-dlon,latr-dlat,datrp,undef_value=0.)
  datrp = env + datrp 
  write(*,'("  vortex TC center clon,clat=           ",f6.1,f6.1)') clon,clat
  write(*,'("  relocated TC center clon_new,clat_new=",f6.1,f6.1)') clon_new,clat_new
  write(*,'("  difference of TC center dlon,dlat=    ",f6.2,f6.2)') dlon,dlat
  write(*,*)
  ! for debug
  !call write_grads(21,'check_plev_relocated.bin',irmax,jrmax,kmax,1,8,lonr,latr,datrp)


  !!
  !! 7) horizontal interpolation to gaussain grid
  !!
  write(*,*) "[7] HORIZONTAL INTERPOLATION FROM REGIONAL GRID TO GAUSSIAN GRID"
  write(*,*)
  !!
  !! need two rows & columns to use quasi-cubic interpolation
  !!
  !                                       ! R(2,jrmax-1)     -        R(irmax-1,jrmax-1)
  ix1 = searchidx(lon,lonr(2),-1)         !       G(ix1,jy2) - G(ix2,jy2)
  ix2 = searchidx(lon,lonr(irmax-1),1)    ! |          |     |              |
  jy1 = searchidx(lat,latr(2),-1)         !       G(ix1,jy1) - G(ix2,jy1)
  jy2 = searchidx(lat,latr(jrmax-1),1)    ! R(2,2)           -        R(irmax-1,2)
  irgmax = ix2-ix1+1                      !   
  jrgmax = jy2-jy1+1                      !
  ijrgmax = irgmax*jrgmax

  deallocate( datrgp )
  allocate( datrgp(irgmax,jrgmax,nvmax) )
  call interp2d_quasicubic(irmax,jrmax,nvmax,lonr,latr,datrp,irgmax,jrgmax,lon(ix1:ix2),lat(jy1:jy2),datrgp) 


  !!
  !! 8) vertical interpolation to model level
  !!
  write(*,*) "[8] VERTICAL INTERPOLATION TO THE MODEL LEVEL"  
  write(*,*)
  deallocate( prgm, prgp )
  allocate( prgm(irgmax,jrgmax,kmax), prgp(irgmax,jrgmax,kmax ) )

  !! 8.1) calculate surface pressure
  do k = 1, kmax
     prgp(:,:,k) = plev(k)
  end do
  psg(ix1:ix2,jy1:jy2) = calcmet_slp2ps(datrgp(:,:,nvt:nvt+dk),hsg(ix1:ix2,jy1:jy2),datrgp(:,:,nvslp),prgp)
  !! 8.2) calculate model pressure level corespoding to the relocated field
  call sigio_modpr(ijrgmax,ijrgmax,kmax,nvcoord,idvc,idsl,vcoord,ierr,ps=psg(ix1:ix2,jy1:jy2),pm=prgm)

  !! 8.3) vertical interpolation
  nv = nvt
  call p2p_extrapolate_T(ijrgmax,kmax,datrgp(:,:,nv:nv+dk),prgp,psg(ix1:ix2,jy1:jy2),  &
       &                 kmax,prgm,psg(ix1:ix2,jy1:jy2),hsg(ix1:ix2,jy1:jy2),tg(ix1:ix2,jy1:jy2,:))
  nv = nv + kmax
  call p2p(ijrgmax,kmax,datrgp(:,:,nv:nv+dk),prgp,psg(ix1:ix2,jy1:jy2),kmax,prgm,dg(ix1:ix2,jy1:jy2,:)) 
  nv = nv + kmax
  call p2p(ijrgmax,kmax,datrgp(:,:,nv:nv+dk),prgp,psg(ix1:ix2,jy1:jy2),kmax,prgm,zg(ix1:ix2,jy1:jy2,:))
  nv = nvq
  call p2p_q2rh(ijrgmax,kmax,datrgp(:,:,nv:nv+dk),datrgp(:,:,nvt:nvt+dk),prgp,psg(ix1:ix2,jy1:jy2), &
       &        kmax,prgm,psg(ix1:ix2,jy1:jy2),hsg(ix1:ix2,jy1:jy2),qg(ix1:ix2,jy1:jy2,:)) 
  nv = nv + kmax
  call p2p(ijrgmax,kmax,datrgp(:,:,nv:nv+dk),prgp,psg(ix1:ix2,jy1:jy2),kmax,prgm,o3g(ix1:ix2,jy1:jy2,:)) 
  nv = nv + kmax
  call p2p(ijrgmax,kmax,datrgp(:,:,nv:nv+dk),prgp,psg(ix1:ix2,jy1:jy2),kmax,prgm,cwcg(ix1:ix2,jy1:jy2,:)) 


  !!
  !! 9) write output gfs sigma file
  !!
  write(*,*) "[9] WRITE OUTPUT GFS SIGMA FILE"
  call gfs_grid2sigma(imax,jmax,kmax,maxwv,4,headi,datao,hsg,psg,tg,dg,zg,qg,o3g,cwcg,yrev=.false.)
  call gfs_sigmawrite(fno,'output.sigma',headi,datao)


end program main




