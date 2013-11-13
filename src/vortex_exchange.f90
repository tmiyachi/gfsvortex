program main
  use type_module
  use constant_module, only : PI=>math_pi  
  use sigio_module
  use gfs_module
  use p2p_module
  use calcmet_module
  use ip_module
  use vortex_module
  use debug_module
  implicit none
  !!
  ! PROGRAM vortex_exchange
  !
  ! DESCRIPTION
  !  separate two input GFS sigma files into vortex field and environmental field,
  !  merge the vortex field into the another vortex field and write out GFS sigma file.
  !  This program contains following steps
  !   1) interpolate golobal field data into 40x40 degrees subdomain area
  !   2) vertically interpolate the field from model level to constant pressure level
  !   3) separate into vortex field and environmental field
  !   4) merge the vortex field and the other environmental field
  !   5) vertically interpolate the field from constant pressure level to model level
  !   6) interpolate subdomain grid to global grid
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
  !
  ! LIBRARY
  !  sigio - sigio_rrohdc,sigio_cnvtdv,sigio_modpr
  !  splib - splat,sptez,sptezm,sptezmv
  !
  ! AUTHOR
  !  Tetsuro Miyachi
  !
  !!
  ! PARAMETER - irmax,jrmax,dxyr
  !  size of subdomain grid to applay vortex removal scheme
  !  41x41 horizontal grids with grid intervals of 1 degree (Liu et al. (2000))
  integer(kind=i4b), parameter :: irmax=41,jrmax=41,ijrmax=irmax*jrmax
  real(kind=sp),     parameter :: dxyr=1.

  integer(kind=i4b), parameter :: fnm=11,fni1=12,fni2=13,fno=21

  integer(kind=i4b) :: imax1,jmax1,kmax1
  integer(kind=i4b) :: maxwv1,ntrac1,idsl1,idvc1,idvm1,nvcoord1,idate1(4)
  integer(kind=i4b) :: imax2,jmax2,kmax2
  integer(kind=i4b) :: maxwv2,ntrac2,idsl2,idvc2,idvm2,nvcoord2,idate2(4)
  type(sigio_head) :: headi1,headi2,heado
  type(sigio_data) :: datai1,datai2,datao
  real(kind=sp), allocatable :: vcoord1(:,:),vcoord2(:,:)

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

  real(kind=sp)              :: lonr(irmax),latr(jrmax)
  real(kind=sp), allocatable :: hsr1(:,:),psr1(:,:),tr1(:,:,:)
  real(kind=sp), allocatable :: dr1(:,:,:),zr1(:,:,:),ur1(:,:,:),vr1(:,:,:)
  real(kind=sp), allocatable :: qr1(:,:,:),o3r1(:,:,:),cwcr1(:,:,:)
  real(kind=sp), allocatable :: hsr2(:,:),psr2(:,:),tr2(:,:,:)
  real(kind=sp), allocatable :: dr2(:,:,:),zr2(:,:,:),ur2(:,:,:),vr2(:,:,:)
  real(kind=sp), allocatable :: qr2(:,:,:),o3r2(:,:,:),cwcr2(:,:,:)
  real(kind=sp), allocatable :: prm1(:,:,:),prm2(:,:,:),datrp1(:,:,:),datrp2(:,:,:)
  real(kind=sp), allocatable :: psref(:,:),prp(:,:,:)

  real(kind=sp), allocatable :: lu(:,:),lv(:,:)
  real(kind=sp), allocatable :: env1(:,:,:),vrtex1(:,:,:),env2(:,:,:),vrtex2(:,:,:)

  integer(kind=i4b) :: klu,klv,ix1,ix2,jy1,jy2,ixc,jyc

  integer(kind=i4b) :: irgmax,jrgmax
  real(kind=sp)     :: clon1,clat1,clon2,clat2,dlon,dlat
  integer(kind=i4b) :: i,j,k,n,ierr,dk,m
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
  !! 2) interpolate to the regional subdomain 
  !!
  write(*,*) "[2] INTERPOLATION TO THE REGIONAL SUBDOMAIN"

  ixc = searchidx(lon2,clon_obs,0)    
  do i = 1, irmax
     lonr(i) = lon2(ixc) + dxyr*real(-irmax/2 + (i-1))
  end do
  jyc = searchidx(lat2,clat_obs,0)    
  do j = 1, jrmax
     latr(j) = lat2(jyc) + dxyr*real(-jrmax/2 + (j-1))
  end do
  write(*,'(i4" x ",i0," domain ",f4.1," degrees interval")'),irmax,jrmax,dxyr
  write(*,'("  domain center is     ",f7.2,f7.2)'),lon2(ixc),lat2(jyc)
  write(*,'("  lonr(1),lonr(irmax)= ",f7.2,f7.2)'),lonr(1),lonr(irmax)
  write(*,'("  latr(1),latr(jrmax)= ",f7.2,f7.2)'),latr(1),latr(jrmax)
  write(*,*)

  allocate( hsr1(irmax,jrmax), psr1(irmax,jrmax))
  allocate( tr1(irmax,jrmax,kmax1) )
  allocate( dr1(irmax,jrmax,kmax1), zr1(irmax,jrmax,kmax1) )
  allocate( ur1(irmax,jrmax,kmax1), vr1(irmax,jrmax,kmax1) )
  allocate( qr1(irmax,jrmax,kmax1), o3r1(irmax,jrmax,kmax1), cwcr1(irmax,jrmax,kmax1) )
  allocate( hsr2(irmax,jrmax), psr2(irmax,jrmax))
  allocate( tr2(irmax,jrmax,kmax2) )
  allocate( dr2(irmax,jrmax,kmax2), zr2(irmax,jrmax,kmax2) )
  allocate( ur2(irmax,jrmax,kmax2), vr2(irmax,jrmax,kmax2) )
  allocate( qr2(irmax,jrmax,kmax2), o3r2(irmax,jrmax,kmax2), cwcr2(irmax,jrmax,kmax2) )

  !! 2.1) environmental field
  call interp2d_quasicubic(imax1,jmax1,1,lon1,lat1,hsg1,irmax,jrmax,lonr,latr,hsr1)
  call interp2d_quasicubic(imax1,jmax1,1,lon1,lat1,psg1,irmax,jrmax,lonr,latr,psr1)
  call interp2d_quasicubic(imax1,jmax1,kmax1,lon1,lat1,tg1,irmax,jrmax,lonr,latr,tr1)
  call interp2d_quasicubic(imax1,jmax1,kmax1,lon1,lat1,dg1,irmax,jrmax,lonr,latr,dr1)
  call interp2d_quasicubic(imax1,jmax1,kmax1,lon1,lat1,zg1,irmax,jrmax,lonr,latr,zr1)
  call interp2d_quasicubic(imax1,jmax1,kmax1,lon1,lat1,ug1,irmax,jrmax,lonr,latr,ur1)
  call interp2d_quasicubic(imax1,jmax1,kmax1,lon1,lat1,vg1,irmax,jrmax,lonr,latr,vr1)
  call interp2d_quasicubic(imax1,jmax1,kmax1,lon1,lat1,qg1,irmax,jrmax,lonr,latr,qr1)
  call interp2d_quasicubic(imax1,jmax1,kmax1,lon1,lat1,o3g1,irmax,jrmax,lonr,latr,o3r1)
  call interp2d_quasicubic(imax1,jmax1,kmax1,lon1,lat1,cwcg1,irmax,jrmax,lonr,latr,cwcr1)

  !! 2.2) vortex field
  call interp2d_quasicubic(imax2,jmax2,2,lon2,lat2,hsg2,irmax,jrmax,lonr,latr,hsr2)
  call interp2d_quasicubic(imax2,jmax2,2,lon2,lat2,psg2,irmax,jrmax,lonr,latr,psr2)
  call interp2d_quasicubic(imax2,jmax2,kmax2,lon2,lat2,tg2,irmax,jrmax,lonr,latr,tr2)
  call interp2d_quasicubic(imax2,jmax2,kmax2,lon2,lat2,dg2,irmax,jrmax,lonr,latr,dr2)
  call interp2d_quasicubic(imax2,jmax2,kmax2,lon2,lat2,zg2,irmax,jrmax,lonr,latr,zr2)
  call interp2d_quasicubic(imax2,jmax2,kmax2,lon2,lat2,ug2,irmax,jrmax,lonr,latr,ur2)
  call interp2d_quasicubic(imax2,jmax2,kmax2,lon2,lat2,vg2,irmax,jrmax,lonr,latr,vr2)
  call interp2d_quasicubic(imax2,jmax2,kmax2,lon2,lat2,qg2,irmax,jrmax,lonr,latr,qr2)
  call interp2d_quasicubic(imax2,jmax2,kmax2,lon2,lat2,o3g2,irmax,jrmax,lonr,latr,o3r2)
  call interp2d_quasicubic(imax2,jmax2,kmax2,lon2,lat2,cwcg2,irmax,jrmax,lonr,latr,cwcr2)

  !for debug
  !call write_grads(21,'check_psr1_total.bin',irmax,jrmax,kmax1,1,0,lonr,latr,psr1)  
  !call write_grads(21,'check_psr2_total.bin',irmax,jrmax,kmax2,1,0,lonr,latr,psr2)  



  !!
  !! 3) vertical interpolate to the constant pressure level
  !!
  write(*,*) "[3] VERTICAL INTERPOLATION TO THE CONSTANT PRESSURE LEVEL"

  !! 3.1) calculate input model pressure level
  allocate( prm1(irmax,jrmax,kmax1), prm2(irmax,jrmax,kmax2) )
  allocate( psref(irmax,jrmax), prp(irmax,jrmax,kmax1) )
  call sigio_modpr(ijrmax,ijrmax,kmax1,nvcoord1,idvc1,idsl1,vcoord1,ierr,ps=psr1,pm=prm1)
  call sigio_modpr(ijrmax,ijrmax,kmax2,nvcoord2,idvc2,idsl2,vcoord2,ierr,ps=psr2,pm=prm2)

  !! 3.2) calculate constant pressure level for vertical interpolation
  psref = max( maxval(psr1), maxval(psr2) )    ! referrence pressure
  call sigio_modpr(ijrmax,ijrmax,kmax1,nvcoord1,idvc1,idsl1,vcoord1,ierr,ps=psref,pm=prp)

  !! 3.3) vertically interpolate data from model levels to constant pressure levels
  nvmax = 1 + 5*kmax1 + 3*kmax1
  dk = kmax1 - 1
  allocate( datrp1(irmax,jrmax,nvmax), datrp2(irmax,jrmax,nvmax) )

  !! environmental field
  nv = 2
  nvt = nv
  call p2p_extrapolate_T(ijrmax,kmax1,tr1,prm1,psr1,kmax1,prp,psr1,hsr1,datrp1(:,:,nv:nv+dk))
  nv = nv + kmax1
  call p2p(ijrmax,kmax1,dr1,prm1,psr1,kmax1,prp,datrp1(:,:,nv:nv+dk))
  nv = nv + kmax1
  call p2p(ijrmax,kmax1,zr1,prm1,psr1,kmax1,prp,datrp1(:,:,nv:nv+dk))
  nv = nv + kmax1
  nvu = nv
  call p2p(ijrmax,kmax1,ur1,prm1,psr1,kmax1,prp,datrp1(:,:,nv:nv+dk))
  nv = nv + kmax1
  nvv = nv
  call p2p(ijrmax,kmax1,vr1,prm1,psr1,kmax1,prp,datrp1(:,:,nv:nv+dk))
  nv = nv + kmax1
  nvq = nv
  qr1 = calcmet_rh(tr1,qr1,prm1)   !SPH->RH        
  call p2p(ijrmax,kmax1,qr1,prm1,psr1,kmax1,prp,datrp1(:,:,nv:nv+dk))
  datrp1(:,:,nv:nv+dk) = max(min(datrp1(:,:,nv:nv+dk),100.),0.)
  datrp1(:,:,nv:nv+dk) = calcmet_q(datrp1(:,:,nvt:nvt+dk),datrp1(:,:,nv:nv+dk),prp)    !RH->SPH
  nv = nv + kmax1
  call p2p(ijrmax,kmax1,o3r1,prm1,psr1,kmax1,prp,datrp1(:,:,nv:nv+dk))
  nv = nv + kmax1
  call p2p(ijrmax,kmax1,cwcr1,prm1,psr1,kmax1,prp,datrp1(:,:,nv:nv+dk))
    
  ! vortex field
  nv = nvt
  call p2p_extrapolate_T(ijrmax,kmax2,tr2,prm2,psr2,kmax1,prp,psr2,hsr2,datrp2(:,:,nv:nv+dk))
  nv = nv + kmax1
  call p2p(ijrmax,kmax2,dr2,prm2,psr2,kmax1,prp,datrp2(:,:,nv:nv+dk))
  nv = nv + kmax1
  call p2p(ijrmax,kmax2,zr2,prm2,psr2,kmax1,prp,datrp2(:,:,nv:nv+dk))
  nv = nv + kmax1
  call p2p(ijrmax,kmax2,ur2,prm2,psr2,kmax1,prp,datrp2(:,:,nv:nv+dk))
  nv = nv + kmax1
  call p2p(ijrmax,kmax2,vr2,prm2,psr2,kmax1,prp,datrp2(:,:,nv:nv+dk))
  nv = nv + kmax1
  qr2 = calcmet_rh(tr2,qr2,prm2)   !SPH->RH        
  call p2p(ijrmax,kmax2,qr2,prm2,psr2,kmax1,prp,datrp2(:,:,nv:nv+dk))
  datrp2(:,:,nv:nv+dk) = max(min(datrp2(:,:,nv:nv+dk),100.),0.)
  datrp2(:,:,nv:nv+dk) = calcmet_q(datrp2(:,:,nvt:nvt+dk),datrp2(:,:,nv:nv+dk),prp)    !RH->SPH
  nv = nv + kmax1
  call p2p(ijrmax,kmax2,o3r2,prm2,psr2,kmax1,prp,datrp2(:,:,nv:nv+dk))
  nv = nv + kmax1
  call p2p(ijrmax,kmax2,cwcr2,prm2,psr2,kmax1,prp,datrp2(:,:,nv:nv+dk))

  !! 3.4) convert surface pressure to sea level pressure
  nvslp = 1
  datrp1(:,:,nvslp) = calcmet_slp(datrp1(:,:,nvt:nvt+dk),hsr1,psr1,prp) 
  datrp2(:,:,nvslp) = calcmet_slp(datrp2(:,:,nvt:nvt+dk),hsr2,psr2,prp) 

  write(*,'("  refference pressure is ",f9.2)') psref(1,1) 
  write(*,'("  prm1(1,1,1)= ",f9.2," prm2(1,1,1)= ",f9.2," prp(1,1,1)= ",f9.2)') prm1(1,1,1),prm2(1,1,1),prp(1,1,1)
  write(*,*)


  !!
  !! 4) separate vortex field
  !!
  write(*,*) "[4] COMPUTE VORTEX FIELD AND ENVIRONMENTAL FIELD"
  k = searchidx(prp(1,1,:),85000.,0)
  allocate( lu(irmax,jrmax), lv(irmax,jrmax) )
  allocate( env1(irmax,jrmax,nvmax),vrtex1(irmax,jrmax,nvmax))
  allocate( env2(irmax,jrmax,nvmax),vrtex2(irmax,jrmax,nvmax))
  write(*,'("  determine filter domain at ",f9.2,"Pa, k=",i0)'),prp(1,1,k),k 
  write(*,*)
  lu = datrp1(:,:,nvu+k-1)
  lv = datrp1(:,:,nvv+k-1)

  write(*,*) "[4.1] ENVIRONMENTAL DATA"
  call separate_env_vortex(irmax,jrmax,nvmax,lonr,latr,datrp1,lu,lv,clon_obs,clat_obs,env1,vrtex1,clon1,clat1)
  !for debug
  call separate_env_vortex(irmax,jrmax,nvmax,lonr,latr,datrp1,lu,lv,clon_obs,clat_obs,env1,vrtex1,clon1,clat1,r0domain=r0domain) 
  call write_domain(21,'check_domain1.txt',clon1,clat1,r0domain)
  lu = datrp2(:,:,nvu+k-1)
  lv = datrp2(:,:,nvv+k-1)

  write(*,*) "[4.1] VORTEX DATA"
  !call separate_env_vortex(irmax,jrmax,nvmax,lonr,latr,datrp2,lu,lv,clon_obs,clat_obs,env2,vrtex2,clon2,clat2)
  !for debug
  call separate_env_vortex(irmax,jrmax,nvmax,lonr,latr,datrp2,lu,lv,clon_obs,clat_obs,env2,vrtex2,clon2,clat2,r0domain=r0domain)
  call write_domain(21,'check_domain2.txt',clon2,clat2,r0domain)
  write(*,*)

  !for debug
  !call write_grads(31,'check_plev_total1.bin',irmax,jrmax,kmax1,1,8,lonr,latr,datrp1)
  !call write_grads(31,'check_plev_env1.bin',irmax,jrmax,kmax1,1,8,lonr,latr,env1)
  !call write_grads(31,'check_plev_vrtex1.bin',irmax,jrmax,kmax1,1,8,lonr,latr,vrtex1)
  !call write_grads(31,'check_plev_total2.bin',irmax,jrmax,kmax1,1,8,lonr,latr,datrp2)
  !call write_grads(31,'check_plev_env2.bin',irmax,jrmax,kmax1,1,8,lonr,latr,env2)
  !call write_grads(31,'check_plev_vrtex2.bin',irmax,jrmax,kmax1,1,8,lonr,latr,vrtex2)


 !!
 !! 5) vortex merge
 !!
  write(*,*) "[5] VORTEX MERGE"
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
  !! 6) vertical interpolation to model level
  !!
  write(*,*) "[6] VERTICAL INTERPOLATION TO THE MODEL LEVEL"  
  write(*,*)

  !! 6.1) calculate surface pressure
  psr1 = calcmet_ps_fromslp(datrp1(:,:,nvt:nvt+dk),hsr1,datrp1(:,:,nvslp),prp)

  !! 6.2) calculate model pressure level corespoding to the relocated field
  call sigio_modpr(ijrmax,ijrmax,kmax1,nvcoord1,idvc1,idsl1,vcoord1,ierr,ps=psr1,pm=prm1)

  !! 6.3) vertical interpolation
  nv = nvt
  call p2p_extrapolate_T(ijrmax,kmax1,datrp1(:,:,nv:nv+dk),prp,psr1,kmax1,prm1,psr1,hsr1,tr1)
  nv = nv + kmax1
  call p2p(ijrmax,kmax1,datrp1(:,:,nv:nv+dk),prp,psr1,kmax1,prm1,dr1)
  nv = nv + kmax1
  call p2p(ijrmax,kmax1,datrp1(:,:,nv:nv+dk),prp,psr1,kmax1,prm1,zr1)
  nv = nvq
  datrp1(:,:,nv:nv+dk) = calcmet_rh(datrp1(:,:,nvt:nvt+dk),datrp1(:,:,nv:nv+dk),prp)   !SPH->RH        
  call p2p(ijrmax,kmax1,datrp1(:,:,nv:nv+dk),prp,psr1,kmax1,prm1,qr1)
  qr1 = calcmet_q(tr1,qr1,prm1)    !RH->SPH
  qr1 = max(min(qr1,100.),0.)
  nv = nv + kmax1
  call p2p(ijrmax,kmax1,datrp1(:,:,nv:nv+dk),prp,psr1,kmax1,prm1,o3r1)
  nv = nv + kmax1
  call p2p(ijrmax,kmax1,datrp1(:,:,nv:nv+dk),prp,psr1,kmax1,prm1,cwcr1)
  
  !for debug
  !call write_grads(21,'check_psr_merged.bin',irmax,jrmax,kmax1,1,0,lonr,latr,psr1)  


  !!
  !! 7) interpolation to global grid
  !!
  write(*,*) "[7] INTERPOLATION TO GLOBAL GRID"
  !                                       ! R(2,jrmax-1)      -            R(irmax-1,jrmax-1)
  ix1 = searchidx(lon1,lonr(2),-1)        !        G(ix1,jy2) - G(ix2,iy2)
  ix2 = searchidx(lon1,lonr(irmax-1),1)   !    |        |           |         |
  jy1 = searchidx(lat1,latr(2),-1)        !        G(ix1,jy1) - G(ix2,jy1)
  jy2 = searchidx(lat1,latr(jrmax-1),1)   ! R(2,2)            -            R(irmax-1,2)
  irgmax = ix2-ix1+1                      !   
  jrgmax = jy2-jy1+1                      !

  write(*,'("  lon1(ix1),lon1(ix2),lat1(jy1),lat1(jy2)= ",f7.2,f7.2,f6.2,f6.2)'),lon1(ix1),lon1(ix2),lat1(jy1),lat1(jy2)
  write(*,*) 
  nv = 1  
  call interp2d_quasicubic(irmax,jrmax,1,lonr,latr,psr1,irgmax,jrgmax,lon1(ix1:ix2),lat1(jy1:jy2),psg1(ix1:ix2,jy1:jy2))  
  call interp2d_quasicubic(irmax,jrmax,kmax1,lonr,latr,tr1,irgmax,jrgmax,lon1(ix1:ix2),lat1(jy1:jy2),tg1(ix1:ix2,jy1:jy2,:))  
  call interp2d_quasicubic(irmax,jrmax,kmax1,lonr,latr,dr1,irgmax,jrgmax,lon1(ix1:ix2),lat1(jy1:jy2),dg1(ix1:ix2,jy1:jy2,:))  
  call interp2d_quasicubic(irmax,jrmax,kmax1,lonr,latr,zr1,irgmax,jrgmax,lon1(ix1:ix2),lat1(jy1:jy2),zg1(ix1:ix2,jy1:jy2,:))  
  call interp2d_quasicubic(irmax,jrmax,kmax1,lonr,latr,qr1,irgmax,jrgmax,lon1(ix1:ix2),lat1(jy1:jy2),qg1(ix1:ix2,jy1:jy2,:))  
  if (o3merge) then
     call interp2d_quasicubic(irmax,jrmax,kmax1,lonr,latr,o3r1,irgmax,jrgmax,lon1(ix1:ix2),lat1(jy1:jy2),o3g1(ix1:ix2,jy1:jy2,:))  
  end if
  if (cwcmerge) then
     call interp2d_quasicubic(irmax,jrmax,kmax1,lonr,latr,cwcr1,irgmax,jrgmax,lon1(ix1:ix2),lat1(jy1:jy2),cwcg1(ix1:ix2,jy1:jy2,:))  
  end if

  !!
  !! 8) write output gfs sigma file
  !!
  write(*,*) "WRITE OUTPUT GFS SIGMA FILE"
  call gfs_grid2sigma(imax1,jmax1,kmax1,maxwv1,4,headi1,datao,hsg1,psg1,tg1,dg1,zg1,qg1,o3g1,cwcg1,yrev=.false.)
  call gfs_sigmawrite(fno,'output.sigma',headi1,datao)

end program main




