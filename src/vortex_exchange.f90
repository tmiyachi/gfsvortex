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
  ! PROGRAM vortex_exchange
  !
  ! DESCRIPTION
  !  Tis program separates two input GFS sigma files into vortex field and environmental field,
  !  merges the vortex field into the other environmental field and write out GFS sigma file.
  !  This program contains following steps.
  !
  !   1. convert surface pressure to sea level pressure to avoid the steep montain effects.
  !   2. vertically interpolate from model level to constant pressure level
  !      within regional subdomain.
  !   3. horizontally interpolate the regional field from the gaussian grid to the equally
  !      spaced lat-lon grid on the each constant pressure level.
  !   4. split the total field into environmental field and vortex field on the each constant
  !      pressure level using the scheme proposed by Kurihara et al.(1992,1995)
  !   5. add the vortex field of input file to the environmental field ot the other input file.
  !      relocate the vortex component to the specified position if needed.
  !   6. horizontally interpolate the the regional field from the equally spaced lat-lon grid
  !      to the gaussian grid on the each constant pressure level.
  !   7. convert sea level pressure to surface pressure.
  !   8. vertically interpolate from the constant pressure level to the model level coressponding 
  !      to the relocated surface pressure.
  !
  ! FILES
  !  INPUT:
  !    unit 11: namelist file
  !    unit 12: input gfs sigma file for environmental component
  !    unit 13: input gfs sigma file for vortex component
  !  OUTPUT:
  !    unit 21: output gfs sigma file
  !
  ! NAMELIST
  !  tcinfo
  !    clon_obs, clat_obs : observed TC center position for first guess position
  !  param
  !    clon_new, clat_new : new TC center position of relocated vortex
  !    relocate           : 'obs' - relocate vortex center to (clon_obs,clat_obs)
  !                         'env' - relocate vortex center to vortex center of environmental component
  !                         'new' - relocate vortex center user specified position (clon_new,clat_new)
  !                         'vrt' - no relocation is done 
  !    o3merge,cwcmerge   : flag for wheter ozone and cloud water field are exchanged. 
  !                         If flag is False, environmental components are used within vortex domain.
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
  !  2013.11.13 - initial version
  !  see the commit log of the repositories
  !
  !!
  ! PARAMETER - irmax,jrmax,dxyr
  !  size of subdomain grid to applay vortex removal scheme
  !  41x41 horizontal grids with grid intervals of 1 degree (Liu et al. (2000))
  !!
  integer(kind=i4b), parameter :: fnm=11,fni1=12,fni2=13,fno=21

  integer(kind=i4b), parameter :: irmax=41,jrmax=41
  real(kind=sp),     parameter :: dxyr=1.

  integer(kind=i4b)          :: imax1,jmax1,kmax1
  integer(kind=i4b)          :: maxwv1,ntrac1,idsl1,idvc1,idvm1,nvcoord1,idate1(4)
  integer(kind=i4b)          :: imax2,jmax2,kmax2
  integer(kind=i4b)          :: maxwv2,ntrac2,idsl2,idvc2,idvm2,nvcoord2,idate2(4)
  real(kind=sp), allocatable :: vcoord1(:,:),vcoord2(:,:)
  type(sigio_head)           :: headi1,headi2,heado
  type(sigio_data)           :: datai1,datai2,datao

  real(kind=sp)    :: clon_obs,clat_obs,clon_new,clat_new
  character(len=3) :: relocate
  logical          :: o3merge,cwcmerge

  real(kind=sp), allocatable :: lon1(:),lat1(:),lon2(:),lat2(:)
  real(kind=sp), allocatable :: hsg1(:,:),psg1(:,:),tg1(:,:,:)
  real(kind=sp), allocatable :: dg1(:,:,:),zg1(:,:,:),ug1(:,:,:),vg1(:,:,:)
  real(kind=sp), allocatable :: qg1(:,:,:),o3g1(:,:,:),cwcg1(:,:,:)
  real(kind=sp), allocatable :: hsg2(:,:),psg2(:,:),tg2(:,:,:)
  real(kind=sp), allocatable :: dg2(:,:,:),zg2(:,:,:),ug2(:,:,:),vg2(:,:,:)
  real(kind=sp), allocatable :: qg2(:,:,:),o3g2(:,:,:),cwcg2(:,:,:)

  real(kind=sp)              :: psref(1)
  real(kind=sp), allocatable :: plev(:)

  integer(kind=i4b)          :: irgmax1,jrgmax1,ijrgmax1,irgmax2,jrgmax2,ijrgmax2
  real(kind=sp), allocatable :: datrgp1(:,:,:), prgm1(:,:,:), prgp1(:,:,:)
  real(kind=sp), allocatable :: datrgp2(:,:,:), prgm2(:,:,:), prgp2(:,:,:)

  real(kind=sp)              :: lonr(irmax),latr(jrmax)
  real(kind=sp), allocatable :: datrp1(:,:,:),datrp2(:,:,:)
  real(kind=sp), allocatable :: lu(:,:),lv(:,:)
  real(kind=sp), allocatable :: env1(:,:,:),vrtex1(:,:,:),env2(:,:,:),vrtex2(:,:,:)

  integer(kind=i4b) :: klu,klv,ix11,ix21,jy11,jy21,ix12,ix22,jy12,jy22,ixc,jyc

  real(kind=sp)     :: clon1,clat1,clon2,clat2,dlon,dlat
  integer(kind=i4b) :: i,j,k,n,ierr,dk
  integer(kind=i4b) :: nv,nvmax,nvt,nvu,nvv,nvq,nvslp

  !for debug
  real(kind=sp) :: r0domain(2,24)

  namelist/tcinfo/clon_obs,clat_obs
  namelist/param/clon_new,clat_new,relocate,o3merge,cwcmerge

  !!
  !! 0) read namelist
  !!  
  write(*,*) "[0] READ NAMELIST"
  open(fnm,file="vortex_move_namelist",iostat=ierr)
  if ( ierr /= 0 ) then
     print*, "namelist file open error. ierr=", ierr
     call abort
  end if
  read(fnm,nml=tcinfo)
  clon_new = clon_obs
  clat_new = clat_obs
  relocate = 'env'
  o3merge  = .true.
  cwcmerge = .true.
  read(fnm,nml=param)
  close(fnm)
  write(*,'("  observed TC center clon,clat=  ",f5.1,f5.1)') clon_obs,clat_obs
  write(*,*)


  !!
  !! 1) read global grid data
  !!
  write(*,*) "[1] OPEN INPUT GFS SIGMA FILE"

  !! 1.1) environmental field
  call gfs_sigmaopen(fni1,'inputenv.sigma',headi1,datai1)  
  imax1    = headi1%lonb
  jmax1    = headi1%latb
  kmax1    = headi1%levs
  ntrac1   = headi1%ntrac
  maxwv1   = headi1%jcap
  idvc1    = headi1%idvc
  idsl1    = headi1%idsl
  idvm1    = headi1%idvm
  nvcoord1 = headi1%nvcoord
  idate1   = headi1%idate
  allocate( vcoord1(kmax1+1,nvcoord1))
  vcoord1  = headi1%vcoord

  allocate( lon1(imax1), lat1(jmax1) )
  allocate( hsg1(imax1,jmax1), psg1(imax1,jmax1) )
  allocate( tg1(imax1,jmax1,kmax1) )
  allocate( dg1(imax1,jmax1,kmax1), zg1(imax1,jmax1,kmax1) )
  allocate( ug1(imax1,jmax1,kmax1), vg1(imax1,jmax1,kmax1) )
  allocate( qg1(imax1,jmax1,kmax1), o3g1(imax1,jmax1,kmax1), cwcg1(imax1,jmax1,kmax1) )

  call gfs_sigma2grid(imax1,jmax1,kmax1,maxwv1,4,headi1,datai1,lon1,lat1,hsg1,psg1,tg1, &
       &              dg1,zg1,ug1,vg1,qg1,o3g1,cwcg1,yrev=.false.)

  write(*,*) "  environmental field file"
  write(*,'("    imax,jmax,kmax= ",i5,i5,i5)') imax1,jmax1,kmax1
  write(*,'("    maxwv= ",i0, "  idate=",i5,i3,i3,i3)') maxwv1,idate1(4),idate1(2),idate1(3),idate1(1)

  !! 1.2) vortex field
  call gfs_sigmaopen(fni2,'inputvortex.sigma',headi2,datai2)
  imax2    = headi2%lonb
  jmax2    = headi2%latb
  kmax2    = headi2%levs
  ntrac2   = headi2%ntrac
  maxwv2   = headi2%jcap
  idvc2    = headi2%idvc
  idsl2    = headi2%idsl
  idvm2    = headi2%idvm
  nvcoord2 = headi2%nvcoord
  idate2   = headi2%idate
  allocate( vcoord2(kmax2+1,nvcoord2))
  vcoord2  = headi2%vcoord

  allocate( lon2(imax2), lat2(jmax2) )
  allocate( hsg2(imax2,jmax2), psg2(imax2,jmax2) )
  allocate( tg2(imax2,jmax2,kmax2) )
  allocate( dg2(imax2,jmax2,kmax2), zg2(imax2,jmax2,kmax2) )
  allocate( ug2(imax2,jmax2,kmax2), vg2(imax2,jmax2,kmax2) )
  allocate( qg2(imax2,jmax2,kmax2), o3g2(imax2,jmax2,kmax2), cwcg2(imax2,jmax2,kmax2) )

  call gfs_sigma2grid(imax2,jmax2,kmax2,maxwv2,4,headi2,datai2,lon2,lat2,hsg2,psg2,tg2, &
       &              dg2,zg2,ug2,vg2,qg2,o3g2,cwcg2,yrev=.false.)

  write(*,*) "  vortex field file"
  write(*,'("    imax,jmax,kmax= ",i5,i5,i5)') imax2,jmax2,kmax2
  write(*,'("    maxwv= ",i0, "  idate=",i5,i3,i3,i3)') maxwv2,idate2(4),idate2(2),idate2(3),idate2(1)
  write(*,*) 


  !!
  !! 2) calculate regional subdomain
  !!
  write(*,*) "[2] CALCULATE REGIONAL SUBDOMAIN"
  ixc = searchidx(lon1,clon_obs,0)    
  do i = 1, irmax
     lonr(i) = lon1(ixc) + dxyr*real(-irmax/2 + (i-1))
  end do
  jyc = searchidx(lat1,clat_obs,0)    
  do j = 1, jrmax
     latr(j) = lat1(jyc) + dxyr*real(-jrmax/2 + (j-1))
  end do
  write(*,'(i4" x ",i0," domain ",f4.1," degrees interval")'),irmax,jrmax,dxyr
  write(*,'("  domain center is     ",f7.2,f7.2)'),lon1(ixc),lat1(jyc)
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
  !                                         ! G(ix+1,jy2-1)    -         G(ix2-1,jy2-1)
  ix11 = searchidx(lon1,lonr(1),1) - 1      !       R(1,jrmax) - R(irmax,jrmax)
  ix21 = searchidx(lon1,lonr(irmax),-1) + 1 ! |     |          |                      |
  jy11 = searchidx(lat1,latr(1),1) - 1      !       R(1,1)     - R(irmax,1)
  jy21 = searchidx(lat1,latr(jrmax),-1) + 1 ! G(ix+1,jy1+1)    -         G(ix2-1,jy1+1)
  irgmax1 = ix21-ix11+1                         
  jrgmax1 = jy21-jy11+1                      
  ijrgmax1 = irgmax1*jrgmax1

  ix12 = searchidx(lon2,lonr(1),1) - 1    
  ix22 = searchidx(lon2,lonr(irmax),-1) + 1
  jy12 = searchidx(lat2,latr(1),1) - 1     
  jy22 = searchidx(lat2,latr(jrmax),-1) + 1
  irgmax2 = ix22-ix12+1                         
  jrgmax2 = jy22-jy12+1                      
  ijrgmax2 = irgmax2*jrgmax2

  nvmax = 1 + 5*kmax1 + 3*kmax1
  dk = kmax1 - 1

  allocate( datrgp1(irgmax1,jrgmax1,nvmax), datrgp2(irgmax2,jrgmax2,nvmax) )
  allocate( plev(kmax1) )
  allocate( prgm1(irgmax1,jrgmax1,kmax1), prgp1(irgmax1,jrgmax1,kmax1) )
  allocate( prgm2(irgmax2,jrgmax2,kmax2), prgp2(irgmax2,jrgmax2,kmax1) )

  !! 3.1) calculate input model pressure level
  call sigio_modpr(ijrgmax1,ijrgmax1,kmax1,nvcoord1,idvc1,idsl1,vcoord1,ierr,  &
       &           ps=psg1(ix11:ix21,jy11:jy21),pm=prgm1)
  call sigio_modpr(ijrgmax2,ijrgmax2,kmax2,nvcoord2,idvc2,idsl2,vcoord2,ierr,  &
       &           ps=psg2(ix12:ix22,jy12:jy22),pm=prgm2)

  !! 3.2) calculate constant pressure level for vertical interpolation
  psref = max(maxval(psg1(ix11:ix21,jy11:jy21)), maxval(psg2(ix12:ix22,jy12:jy22))) ! referrence pressure
  call sigio_modpr(1,1,kmax1,nvcoord1,idvc1,idsl1,vcoord1,ierr,ps=psref,pm=plev)
  do k = 1, kmax1
     prgp1(:,:,k) = plev(k)
     prgp2(:,:,k) = plev(k)
  end do

  !! 3.3) convert surface pressure to sea level pressure
  nvslp = 1
  datrgp1(:,:,nvslp) = calcmet_ps2slp(tg1(ix11:ix21,jy11:jy21,:),hsg1(ix11:ix21,jy11:jy21),psg1(ix11:ix21,jy11:jy21),prgm1) 
  datrgp2(:,:,nvslp) = calcmet_ps2slp(tg2(ix12:ix22,jy12:jy22,:),hsg2(ix12:ix22,jy12:jy22),psg2(ix12:ix22,jy12:jy22),prgm2) 

  !! 3.4) vertically interpolate data from model levels to constant pressure levels
  !! environmental component
  nv  = 2
  nvt = nv
  call p2p_extrapolate_T(ijrgmax1,kmax1,tg1(ix11:ix21,jy11:jy21,:),prgm1,psg1(ix11:ix21,jy11:jy21),   &
       &            kmax1,prgp1,psg1(ix11:ix21,jy11:jy21),hsg1(ix11:ix21,jy11:jy21),datrgp1(:,:,nv:nv+dk))
  nv = nv + kmax1
  call p2p(ijrgmax1,kmax1,dg1(ix11:ix21,jy11:jy21,:),prgm1,psg1(ix11:ix21,jy11:jy21),kmax1,prgp1,datrgp1(:,:,nv:nv+dk))
  nv = nv + kmax1
  call p2p(ijrgmax1,kmax1,zg1(ix11:ix21,jy11:jy21,:),prgm1,psg1(ix11:ix21,jy11:jy21),kmax1,prgp1,datrgp1(:,:,nv:nv+dk))
  nv = nv + kmax1
  nvu = nv
  call p2p(ijrgmax1,kmax1,ug1(ix11:ix21,jy11:jy21,:),prgm1,psg1(ix11:ix21,jy11:jy21),kmax1,prgp1,datrgp1(:,:,nv:nv+dk))
  nv = nv + kmax1
  nvv = nv
  call p2p(ijrgmax1,kmax1,vg1(ix11:ix21,jy11:jy21,:),prgm1,psg1(ix11:ix21,jy11:jy21),kmax1,prgp1,datrgp1(:,:,nv:nv+dk))
  nv = nv + kmax1
  nvq = nv  
  call p2p_q2rh(ijrgmax1,kmax1,qg1(ix11:ix21,jy11:jy21,:),tg1(ix11:ix21,jy11:jy21,:),prgm1,psg1(ix11:ix21,jy11:jy21),  &
       &        kmax1,prgp1,psg1(ix11:ix21,jy11:jy21),hsg1(ix11:ix21,jy11:jy21),datrgp1(:,:,nv:nv+dk))  
  nv = nv + kmax1
  call p2p(ijrgmax1,kmax1,o3g1(ix11:ix21,jy11:jy21,:),prgm1,psg1(ix11:ix21,jy11:jy21),kmax1,prgp1,datrgp1(:,:,nv:nv+dk))
  nv = nv + kmax1
  call p2p(ijrgmax1,kmax1,cwcg1(ix11:ix21,jy11:jy21,:),prgm1,psg1(ix11:ix21,jy11:jy21),kmax1,prgp1,datrgp1(:,:,nv:nv+dk))

  !! vortex component
  nv = nvt
  call p2p_extrapolate_T(ijrgmax2,kmax2,tg2(ix12:ix22,jy12:jy22,:),prgm2,psg2(ix12:ix22,jy12:jy22),   &
       &                 kmax1,prgp2,psg2(ix12:ix22,jy12:jy22),hsg2(ix12:ix22,jy12:jy22),datrgp2(:,:,nv:nv+dk))
  nv = nv + kmax1
  call p2p(ijrgmax2,kmax2,dg2(ix12:ix22,jy12:jy22,:),prgm2,psg2(ix12:ix22,jy12:jy22),kmax1,prgp2,datrgp2(:,:,nv:nv+dk))
  nv = nv + kmax1
  call p2p(ijrgmax2,kmax2,zg2(ix12:ix22,jy12:jy22,:),prgm2,psg2(ix12:ix22,jy12:jy22),kmax1,prgp2,datrgp2(:,:,nv:nv+dk))
  nv = nv + kmax1
  call p2p(ijrgmax2,kmax2,ug2(ix12:ix22,jy12:jy22,:),prgm2,psg2(ix12:ix22,jy12:jy22),kmax1,prgp2,datrgp2(:,:,nv:nv+dk))
  nv = nv + kmax1
  call p2p(ijrgmax2,kmax2,vg2(ix12:ix22,jy12:jy22,:),prgm2,psg2(ix12:ix22,jy12:jy22),kmax1,prgp2,datrgp2(:,:,nv:nv+dk))
  nv = nv + kmax1
  call p2p_q2rh(ijrgmax2,kmax2,qg2(ix12:ix22,jy12:jy22,:),tg2(ix12:ix22,jy12:jy22,:),prgm2,psg2(ix12:ix22,jy12:jy22),  &
       &        kmax1,prgp2,psg2(ix12:ix22,jy12:jy22),hsg2(ix12:ix22,jy12:jy22),datrgp2(:,:,nv:nv+dk))  
  nv = nv + kmax1
  call p2p(ijrgmax2,kmax2,o3g2(ix12:ix22,jy12:jy22,:),prgm2,psg2(ix12:ix22,jy12:jy22),kmax1,prgp2,datrgp2(:,:,nv:nv+dk))
  nv = nv + kmax1
  call p2p(ijrgmax2,kmax2,cwcg2(ix12:ix22,jy12:jy22,:),prgm2,psg2(ix12:ix22,jy12:jy22),kmax1,prgp2,datrgp2(:,:,nv:nv+dk))

  write(*,'("  refference pressure is ",f9.2)') psref(1) 
  write(*,'("  prgm1(1,1,1)=       ",f9.2,"  prgp1(1,1,1)=     ",f9.2)') prgm1(1,1,1),prgp1(1,1,1)
  write(*,'("  datrgp1(1,1,nvslp)= ",f9.2,"  datrgp1(1,1,nvt)= ",f9.2)') datrgp1(1,1,nvslp),datrgp1(1,1,nvt)
  write(*,'("  datrgp1(1,1,nvu)=   ",f9.2,"  datrgp1(1,1,nvv)= ",f9.2)') datrgp1(1,1,nvu),datrgp1(1,1,nvv)
  write(*,'("  datrgp1(1,1,nvq)=   ",f9.2)') datrgp1(1,1,nvq)
  write(*,'("  prgm2(1,1,1)=       ",f9.2,"  prgp2(1,1,1)=     ",f9.2)') prgm2(1,1,1),prgp2(1,1,1)
  write(*,'("  datrgp2(1,1,nvslp)= ",f9.2,"  datrgp2(1,1,nvt)= ",f9.2)') datrgp2(1,1,nvslp),datrgp2(1,1,nvt)
  write(*,'("  datrgp2(1,1,nvu)=   ",f9.2,"  datrgp2(1,1,nvv)= ",f9.2)') datrgp2(1,1,nvu),datrgp2(1,1,nvv)
  write(*,'("  datrgp2(1,1,nvq)=   ",f9.2)') datrgp2(1,1,nvq)
  write(*,*)


  !!
  !! 4) horizontal interpolation
  !!
  write(*,*) "[4] HORIZONTAL INTERPOLATION FROM GAUSSIAN GRID TO REGIONAL GRID"
  write(*,*) 
  allocate( datrp1(irmax,jrmax,nvmax), datrp2(irmax,jrmax,nvmax) )
  call interp2d_quasicubic(irgmax1,jrgmax1,nvmax,lon1(ix11:ix21),lat1(jy11:jy21),datrgp1,irmax,jrmax,lonr,latr,datrp1)
  call interp2d_quasicubic(irgmax2,jrgmax2,nvmax,lon2(ix12:ix22),lat2(jy12:jy22),datrgp2,irmax,jrmax,lonr,latr,datrp2)


  !!
  !! 5) separate vortex field
  !!
  write(*,*) "[5] COMPUTE VORTEX FIELD AND ENVIRONMENTAL FIELD"
  k = searchidx(plev,85000.,0)
  allocate( lu(irmax,jrmax), lv(irmax,jrmax) )
  allocate( env1(irmax,jrmax,nvmax),vrtex1(irmax,jrmax,nvmax))
  allocate( env2(irmax,jrmax,nvmax),vrtex2(irmax,jrmax,nvmax))
  write(*,'("  determine filter domain at ",f9.2,"Pa, k=",i0)'),prgp1(1,1,k),k 
  write(*,*)

  write(*,*) "[5.1] ENVIRONMENTAL DATA"
  lu = datrp1(:,:,nvu+k-1)
  lv = datrp1(:,:,nvv+k-1)
  call separate_env_vortex(irmax,jrmax,nvmax,lonr,latr,datrp1,lu,lv,clon_obs,clat_obs,env1,vrtex1,clon1,clat1)
  !for debug
  !call separate_env_vortex(irmax,jrmax,nvmax,lonr,latr,datrp1,lu,lv,clon_obs,clat_obs,env1,vrtex1,clon1,clat1,r0domain=r0domain) 
  !call write_domain(21,'check_domain1.txt',clon1,clat1,r0domain)

  write(*,*) "[5.1] VORTEX DATA"
  lu = datrp2(:,:,nvu+k-1)
  lv = datrp2(:,:,nvv+k-1)
  call separate_env_vortex(irmax,jrmax,nvmax,lonr,latr,datrp2,lu,lv,clon_obs,clat_obs,env2,vrtex2,clon2,clat2)
  !for debug
  !call separate_env_vortex(irmax,jrmax,nvmax,lonr,latr,datrp2,lu,lv,clon_obs,clat_obs,env2,vrtex2,clon2,clat2,r0domain=r0domain)
  !call write_domain(21,'check_domain2.txt',clon2,clat2,r0domain)
  write(*,*)

  !for debug
  !call write_grads(31,'check_plev_total1.bin',irmax,jrmax,kmax1,1,8,lonr,latr,datrp1)
  !call write_grads(31,'check_plev_env1.bin',irmax,jrmax,kmax1,1,8,lonr,latr,env1)
  !call write_grads(31,'check_plev_vrtex1.bin',irmax,jrmax,kmax1,1,8,lonr,latr,vrtex1)
  !call write_grads(31,'check_plev_total2.bin',irmax,jrmax,kmax1,1,8,lonr,latr,datrp2)
  !call write_grads(31,'check_plev_env2.bin',irmax,jrmax,kmax1,1,8,lonr,latr,env2)
  !call write_grads(31,'check_plev_vrtex2.bin',irmax,jrmax,kmax1,1,8,lonr,latr,vrtex2)


 !!
 !! 6) vortex merge
 !!
  write(*,*) "[6] VORTEX MERGE"
  if ( relocate == 'env' ) then
     write(*,*) " relocate the vortex center to the environmental data TC center"
     clon_new = clon1
     clat_new = clat1  
     dlon = clon_new - clon2
     dlat = clat_new - clat2
     call interp2d_quasicubic(irmax,jrmax,nvmax,lonr,latr,vrtex2,irmax,jrmax,lonr-dlon,latr-dlat,vrtex1,undef_value=0.)
     datrp1 = env1 + vrtex1
  else if ( relocate == 'obs' ) then
     write(*,*) " relocate the vortex center to the observed TC center"
     clon_new = clon_obs
     clat_new = clat_obs
     dlon = clon_new - clon2
     dlat = clat_new - clat2
     call interp2d_quasicubic(irmax,jrmax,nvmax,lonr,latr,vrtex2,irmax,jrmax,lonr-dlon,latr-dlat,vrtex1,undef_value=0.)
     datrp1 = env1 + vrtex1
  else if ( relocate == 'new' ) then
     write(*,*) " relocate the vortex center to the user specified TC center"
     dlon = clon_new - clon2
     dlat = clat_new - clat2
     call interp2d_quasicubic(irmax,jrmax,nvmax,lonr,latr,vrtex2,irmax,jrmax,lonr-dlon,latr-dlat,vrtex1,undef_value=0.)
     datrp1 = env1 + vrtex1
  else if ( relocate == 'vrt' ) then
     write(*,*) " no relocation is done"
     clon_new = clon2
     clat_new = clat2
     dlon = 0
     dlat = 0
     datrp1 = env1 + vrtex2
  else
     write(*,*) "invalid relocate option", relocate
     call abort
  end if
  write(*,'("  vortex TC center clon,clat=           ",f6.1,f6.1)') clon2,clat2
  write(*,'("  relocated TC center clon_new,clat_new=",f6.1,f6.1)') clon_new,clat_new
  write(*,'("  difference of TC center dlon,dlat=    ",f6.2,f6.2)') dlon,dlat
  write(*,*)
  
  !for debug
  !call write_grads(31,'check_plev_merged.bin',irmax,jrmax,kmax1,1,8,lonr,latr,datrp1)


  !!
  !! 7) horizontal interpolation to gaussain grid
  !!
  write(*,*) "[7] HORIZONTAL INTERPOLATION FROM REGIONAL GRID TO GAUSSIAN GRID"
  write(*,*)
  !!
  !! need two rows & columns to use quasi-cubic interpolation
  !!
  !                                         ! R(2,jrmax-1)     -        R(irmax-1,jrmax-1)
  ix11 = searchidx(lon1,lonr(2),-1)         !       G(ix1,jy2) - G(ix2,jy2)
  ix21 = searchidx(lon1,lonr(irmax-1),1)    ! |          |     |              |
  jy11 = searchidx(lat1,latr(2),-1)         !       G(ix1,jy1) - G(ix2,jy1)
  jy21 = searchidx(lat1,latr(jrmax-1),1)    ! R(2,2)           -        R(irmax-1,2)
  irgmax1 = ix21-ix11+1     
  jrgmax1 = jy21-jy11+1      
  ijrgmax1 = irgmax1*jrgmax1

  deallocate( datrgp1 )
  allocate( datrgp1(irgmax1,jrgmax1,nvmax) )
  call interp2d_quasicubic(irmax,jrmax,nvmax,lonr,latr,datrp1,irgmax1,jrgmax1,lon1(ix11:ix21),lat1(jy11:jy21),datrgp1) 


  !!
  !! 8) vertical interpolation to model level
  !!
  write(*,*) "[8] VERTICAL INTERPOLATION TO THE MODEL LEVEL"  
  write(*,*)
  deallocate( prgm1, prgp1 )
  allocate( prgm1(irgmax1,jrgmax1,kmax1), prgp1(irgmax1,jrgmax1,kmax1 ) )

  !! 8.1) calculate surface pressure
  do k = 1, kmax1
     prgp1(:,:,k) = plev(k)
  end do
  psg1(ix11:ix21,jy11:jy21) = calcmet_slp2ps(datrgp1(:,:,nvt:nvt+dk),hsg1(ix11:ix21,jy11:jy21),datrgp1(:,:,nvslp),prgp1)

  !! 8.2) calculate model pressure level corespoding to the relocated field
  call sigio_modpr(ijrgmax1,ijrgmax1,kmax1,nvcoord1,idvc1,idsl1,vcoord1,ierr,ps=psg1(ix11:ix21,jy11:jy21),pm=prgm1)

  !! 8.3) vertical interpolation
  nv = nvt
  call p2p_extrapolate_T(ijrgmax1,kmax1,datrgp1(:,:,nv:nv+dk),prgp1,psg1(ix11:ix21,jy11:jy21),  &
       &                 kmax1,prgm1,psg1(ix11:ix21,jy11:jy21),hsg1(ix11:ix21,jy11:jy21),tg1(ix11:ix21,jy11:jy21,:))
  nv = nv + kmax1
  call p2p(ijrgmax1,kmax1,datrgp1(:,:,nv:nv+dk),prgp1,psg1(ix11:ix21,jy11:jy21),kmax1,prgm1,dg1(ix11:ix21,jy11:jy21,:)) 
  nv = nv + kmax1
  call p2p(ijrgmax1,kmax1,datrgp1(:,:,nv:nv+dk),prgp1,psg1(ix11:ix21,jy11:jy21),kmax1,prgm1,zg1(ix11:ix21,jy11:jy21,:))
  nv = nvq
  call p2p_q2rh(ijrgmax1,kmax1,datrgp1(:,:,nv:nv+dk),datrgp1(:,:,nvt:nvt+dk),prgp1,psg1(ix11:ix21,jy11:jy21), &
       &        kmax1,prgm1,psg1(ix11:ix21,jy11:jy21),hsg1(ix11:ix21,jy11:jy21),qg1(ix11:ix21,jy11:jy21,:)) 
  nv = nv + kmax1
  if (o3merge) then
     call p2p(ijrgmax1,kmax1,datrgp1(:,:,nv:nv+dk),prgp1,psg1(ix11:ix21,jy11:jy21),kmax1,prgm1,o3g1(ix11:ix21,jy11:jy21,:)) 
  end if
  nv = nv + kmax1
  if (cwcmerge) then
     call p2p(ijrgmax1,kmax1,datrgp1(:,:,nv:nv+dk),prgp1,psg1(ix11:ix21,jy11:jy21),kmax1,prgm1,cwcg1(ix11:ix21,jy11:jy21,:)) 
  end if


  !!
  !! 9) write output gfs sigma file
  !!
  write(*,*) "[9] WRITE OUTPUT GFS SIGMA FILE"
  call gfs_grid2sigma(imax1,jmax1,kmax1,maxwv1,4,headi1,datao,hsg1,psg1,tg1,dg1,zg1,qg1,o3g1,cwcg1,yrev=.false.)
  call gfs_sigmawrite(fno,'output.sigma',headi1,datao)

end program main




