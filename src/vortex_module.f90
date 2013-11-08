module vortex_module
  !
  ! DESCRIPTION
  !  This module contains the vortex separtion routines based on the scheme 
  !  proposed by Kurihara et al. (1992,1995)
  !
  ! REFFERENCE
  !  Kurihara et al., 1992: An initialization scheme of hurricane models 
  !    by vortex specification. Mon. Wea. Rev., 121, 2030-2045.
  !  Kurihara et al., 1995: Improvements in the GFDL hurricane prediction system. 
  !    Mon. Wea. Rev., 123, 2791-2801.
  !  Liu et al., Improvements in hurricane initialization and forecasting at 
  !    NCEP with global and regional (GFDL) models, NCEP office note 472
  !
  !  NOTE
  !   Sobroutine 'NNLS' is used to solve non-negative least square problems.
  !   avaiable at http://hesperia.gsfc.nasa.gov/~schmahl/nnls/index.html
  !  
  !
  use type_module
  use constant_module, only : PI=>math_pi
  use ip_module
  !for debug
  use debug_module
  implicit none

  private
  public :: separate_env_vortex

  !! PARMETER FOR VORTEX SEPARATION
  !
  ! parameter to calculate filtered domain
  ! 0.1 degree radial intervals and 24 points on the circle 
  ! out to a distance of 12 degree (~1200km) (KBR95)
  integer(kind=i4b), parameter, private :: rrmax=120,prmax=24
  real(kind=sp)    , parameter, private :: drr=0.1

  ! parameter of subdomain to calculate storm center
  ! 7x7 grid by 1 degree (KBR95)
  ! 0.2 degree radial intervals and 24 points on te circle 
  ! out to a distance of 6 degree (KBR95)
  integer(kind=i4b), parameter, private :: idmax=7,jdmax=7,rdmax=30,pdmax=24
  real(kind=sp)    , parameter, private :: dxyd=1,drd=0.2

  ! parameter of optimum interpolation mothod
  ! 2.5 degree (~250km) (KBR95)
  real(kind=sp),     parameter, private :: dscale2 = (2.5)**2  !(dscale)**2

  !for debug
  logical, private :: debug = .false.
contains

  subroutine separate_env_vortex(imax,jmax,kmax,xr,yr,datar,lur,lvr,cx,cy,env,vrtex,cx_new,cy_new)
    !
    ! DESCRIPTION
    !
    ! ARGUMENTS
    !  INPUT:
    !    xr(imax)              : input regional grid x-coordinate
    !    yr(jmax)              : input regional grid y-coordinate
    !    datar(imax,jmax,kmax) : input regional grid data
    !    cx,cy                 : first-guess TC center position
    !    lur(imax,jmax)        : low-level zonal wind to calculate TC center
    !    lvr(imax,jmax)        : low-level meridional wind to calculate TC center
    !  OUTPUT:
    !    env(imax,jmax,kmax)   : environmental component
    !    vrtex(imax,jmax,kmax) : vortex component
    !
    ! MODULE
    !  ip_module - interp2d,searchidx
    !
    ! MODULE PARAMETER
    !  integer(kind=i4b) :: pmaxr
    !
    integer(kind=i4b), intent(in)  :: imax,jmax,kmax
    real(kind=sp),     intent(in)  :: xr(imax),yr(jmax),datar(imax,jmax,kmax)
    real(kind=sp),     intent(in)  :: lur(imax,jmax),lvr(imax,jmax)
    real(kind=sp),     intent(in)  :: cx,cy
    real(kind=sp),     intent(out) :: vrtex(imax,jmax,kmax),env(imax,jmax,kmax)
    real(kind=sp),     intent(out) :: cx_new,cy_new

    integer(kind=i4b) :: i,j,k,icx,jcy
    real(kind=sp)     :: r0(prmax),lurd(imax,jmax),lvrd(imax,jmax)
    real(kind=sp)     :: dstrb(imax,jmax,kmax)

    print*, "  START VORTEX SEPARATION"
    ! 1) split subdomain field into basic field and disturbance field by smoothing filter
    print*, "1) APLICATION OF THE SMOOTHING FILTER"   
    call smoothing_filter(imax,jmax,kmax,datar,dstrb)
    call smoothing_filter(imax,jmax,1,lur,lurd)
    call smoothing_filter(imax,jmax,1,lvr,lvrd)

    ! for debug
    !    call write_grads(21,'check_plev_dstrb.bin',imax,jmax,64,1,8,xr,yr,dstrb)
    
    ! 2) calculate new TC center to apply cirindrical filter
    !    use low-level disturbance field
    print*, "2) CALCULATION OF THE NEW STORM CENTER"
    call storm_center(imax,jmax,xr,yr,lurd,lvrd,cx,cy,cx_new,cy_new)

    ! 3) determine the filter domain
    print*, "3) DETERMINATION OF THE FILTER DOMAIN"
    call filter_domain(imax,jmax,xr,yr,lurd,lvrd,cx_new,cy_new,r0)

    ! 4) remove vortex
    print*, "4) CONSTRUCTION OF THE ENVIRONMENTAL FIELD"
    call remove_vortex(imax,jmax,kmax,xr,yr,dstrb,cx_new,cy_new,r0,vrtex)
    env = datar - vrtex

    print*, "  END VORTEX SEPARATION"

  end subroutine separate_env_vortex

  subroutine smoothing_filter(imax,jmax,kmax,total,dstrb)
    !
    ! DESCRIPTION
    !  applay smoothing filter to split the total field into the basic field 
    !  and the disturbance field
    !
    ! ARGUMENTS
    !  INPUT:
    !    total(imax,jmax,kmax) : total scalar field
    !  OUTPUT:
    !    dstrb(imax,jmax,kmax) : disturbance filtered field
    !
    ! REFFERENCE
    !  Kurihara et al., 1992, MWR
    !
    integer(kind=i4b), intent(in)  :: imax,jmax,kmax
    real(kind=sp),     intent(in)  :: total(imax,jmax,kmax)
    real(kind=sp),     intent(out) :: dstrb(imax,jmax,kmax)

    integer(kind=i4b), parameter  :: nmax=11
    integer(kind=i4b) :: i,j,k,n
    real(kind=sp) :: fk(nmax),xtmp(imax,2),ytmp(jmax,2),work(imax,jmax)
    real(kind=sp) :: m(nmax)=(/ 2.,3.,4.,2.,5.,6.,7.,2.,8.,9.,2. /)  
    
    ! filtering parameter 
    do n = 1, nmax
       fk(n) = 0.5/(1-cos(2.*PI/m(n)))   ! eq.(3.5) Kurihara et al. (1992)
    end do
    
    do k = 1, kmax
       work = 0.
       ! do zonal filter
       do j = 1, jmax
          xtmp(:,1) = total(:,j,k)          
          xtmp(:,2) = xtmp(:,1)
          do n = 1, nmax
             do i = 2, imax-1
                xtmp(i,2) = xtmp(i,2) + fk(n)*(xtmp(i-1,1)+xtmp(i+1,1)-2*xtmp(i,1))
             end do
             xtmp(:,1) = xtmp(:,2)
          end do
          work(:,j) = xtmp(:,2)
       end do

       ! do meridional filter
       do i = 1, imax
          ytmp(:,1) = work(i,:)
          ytmp(:,2) = ytmp(:,1)
          do n = 1, nmax
             do j = 2, jmax-1
                ytmp(j,2) = ytmp(j,2) + fk(n)*(ytmp(j-1,1)+ytmp(j+1,1)-2*ytmp(j,1))
             end do
             ytmp(:,1) = ytmp(:,2)
          end do
          dstrb(i,:,k) = total(i,:,k) - ytmp(:,2)
       end do
    end do
    
  end subroutine smoothing_filter

  subroutine storm_center(imax,jmax,xi,yi,u,v,cxi,cyi,cxo,cyo)
    !
    ! DESCRIPTION
    !   For each grid point within subdomain centered at the grid point 
    !  nearest the first-guess position (cxi,cyi), calculate TC center 
    !  position (cxo,cyo) which maximize the azimuthally averaged
    !  tangential component of the disturbance wind for each profile
    !   See more detail in KBR95 section 2.c
    !
    ! ARGUMENTS
    !  INPUT:
    !    cxi,cyi      : first-guess TC center position
    !    xi(imax)     : input x-coordinate
    !    yi(jmax)     : input y-coordinate
    !    u(imax,jmax) : low level disturbance zonal wind (850hPa)
    !    v(imax,jmax) : low level disturbance meridional wind (850hPa)
    ! OUTPUT:
    !    cxo,cyo      : new TC center position
    !
    ! MODULE
    !  ip_module - grd2cyl
    !
    ! MODULE PAREMETER
    !  integer(kind=i4b) :: idmax,jdmax,rdmax,pdmax
    !  real(kind=sp)     :: dxyd,drd
    !
    integer(kind=i4b), intent(in)  :: imax,jmax
    real(kind=sp),     intent(in)  :: cxi,cyi
    real(kind=sp),     intent(in)  :: xi(imax),yi(jmax),u(imax,jmax),v(imax,jmax)
    real(kind=sp),     intent(out) :: cxo,cyo

    integer(kind=i4b) :: i,j,n,m
    real(kind=sp)     :: cx,cy,twd_ave,twd_ave_max, r_max, dphi
    real(kind=sp)     :: r(rdmax),phi(pdmax),uc(rdmax,pdmax),vc(rdmax,pdmax)

    ! cyrindrical coordinate
    do i = 1, rdmax
       r(i) = drd*i
    end do
    dphi = 2.*PI/real(pdmax)
    do i = 1, pdmax
       phi(i) = dphi*(i-1)
    end do

    ! search accurate center position whtin subdomain
    twd_ave_max = -10000.
    do i = 1, idmax
       do j = 1, jdmax
          cx = cxi + dxyd*real(-idmax/2 + (i-1)) !assumed center
          cy = cyi + dxyd*real(-jdmax/2 + (j-1)) !assumed center
          
          ! interpolate cylindical coordinate
          call grd2cyl(imax,jmax,1,xi,yi,cx,cy,u,rdmax,pdmax,r,phi,uc)
          call grd2cyl(imax,jmax,1,xi,yi,cx,cy,v,rdmax,pdmax,r,phi,vc)

          do n = 1, rdmax
             twd_ave = 0.
             ! calculate azimuthally averaged tangential wind
             do m = 1, pdmax
                twd_ave = twd_ave + -uc(n,m)*sin(phi(m)) + vc(n,m)*cos(phi(m))
             end do
             twd_ave = twd_ave/real(pdmax)
             if ( twd_ave > twd_ave_max ) then
                cxo = cx
                cyo = cy
                twd_ave_max = twd_ave
                r_max = r(n)
             end if
          end do
       end do
    end do

    write(*,'("  first-guess TC center  clon,clat=          ",f5.1,f5.1)'),cxi,cyi
    write(*,'("  new TC center  clon_new,clat_new=          ",f5.1,f5.1)'),cxo,cyo
    write(*,'("  maximun averaged tangential wind           ",f5.2)'),twd_ave_max
    write(*,'("  radius of maximum averaged tangential wind ",f5.2)'),r_max

  end subroutine storm_center

  subroutine filter_domain(imax,jmax,xi,yi,u,v,cx,cy,r0)
    !
    ! DESCRIPTION
    !   calculate the extent of the filter domain in each of the directions originating
    !  from the vortex center.
    !   see more detail in KRB95 section 2.c
    !
    ! ARGUMENTS
    !  INPUT:
    !    cx,cy        : TC center position
    !    xi(imax)     : input x-coordinate
    !    yi(jmax)     : input y-coordinate
    !    u(imax,jmax) : low level disturbance zonal wind (850hPa)
    !    v(imax,jmax) : low level disturbance meridional wind (850hPa)
    !  OUTPUT:
    !    r0(prmax)    : firter domain bounding radius at each direction
    !
    ! MODULE
    !  ip_module - grd2cyl, searchidx
    !
    ! MODULE PAREMETER
    !  integer(kind=i4b) :: rrmax,prmax
    !  real(kind=sp)     :: drr
    !
    integer(kind=i4b), intent(in) :: imax,jmax    
    real(kind=sp), intent(in)     :: xi(imax),yi(jmax),u(imax,jmax),v(imax,jmax)
    real(kind=sp), intent(in)     :: cx,cy
    real(kind=sp), intent(out)    :: r0(prmax)

    integer(kind=i4b) :: i, j, n, m, ir_strt,ira,irb
    integer(kind=i4b) :: irf(prmax)
    real(kind=sp) :: r(rrmax), phi(prmax),uc(rrmax,prmax),vc(rrmax,prmax)
    real(kind=sp) :: twd(rrmax,prmax),twd_ave(rrmax),rf(prmax)
    real(kind=sp) :: twd_ave_max,r_max,rf_ave,ra,rb,dphi,e1,e2
    real(kind=sp) :: rs(prmax)

    ! cyrindrical coordinate
    do i = 1, rrmax
       r(i) = drr*i
    end do
    dphi = 2.*PI/real(prmax)
    do i = 1, prmax
       phi(i) = dphi*(i-1)
    end do

    ! convert to cyrindrical coordinate
    call grd2cyl(imax,jmax,1,xi,yi,cx,cy,u,rrmax,prmax,r,phi,uc)
    call grd2cyl(imax,jmax,1,xi,yi,cx,cy,v,rrmax,prmax,r,phi,vc)

    ! calculate tangential wind
    do j = 1, prmax
       e1 = sin(phi(j))
       e2 = cos(phi(j))
       do i = 1, rrmax          
          twd(i,j) = -uc(i,j)*e1 + vc(i,j)*e2
       end do
    end do
    twd_ave = 0.
    do i = 1, rrmax
       do j = 1, prmax
          twd_ave(i) = twd_ave(i) + twd(i,j)
       end do
       twd_ave(i) = twd_ave(i)/real(prmax)
    end do

    !! 1) calclulate starting point to search for the filter radius at each direction
    !!    see more detail in KBR95 Appendix A
    !! 1.1) first starting point is 1.5*Rmax_ave
    i = minval(maxloc(twd_ave))
    twd_ave_max = twd_ave(i)
    r_max = r(i)
    ir_strt = min(searchidx(r,1.5*r_max,-1),rrmax)

    !! 1.2) calculate Rf_ave (algorithm is based on KBR95 section 2.c (see fig.3))
    rf_ave = r_max  ! neither condition is met 
    n = 0
    do i = ir_strt, rrmax
       ! V < 3 m/s
       if ( twd_ave(i) < 3. ) then ! second condition is met
          rf_ave = r(i)
          exit
       end if
       ! V < 6 m/s and -dV/dr < 4x10^-6 /s   !1degree ~ 100km = 1e5 m
       if ( twd_ave(i) < 6. .and. -(twd_ave(i+1)-twd_ave(i))/1.e5/drr < 4.e-6 ) then 
          n = n + 1
          if ( n== 2 ) then       ! first condition is met
             rf_ave = r(i)
             exit
          end if
       end if
    end do
  
    !! 1.3) calculate Rm
    !! parameter: a=0.5,b=0.75 (KBR95)
    ra = 0.5*r_max
    rb = 0.75*r_max + (1-0.75)*rf_ave
    ira = searchidx(r,ra,-1)
    irb = searchidx(r,rb,1)

    !! 2) calculate filtering radius R0=1.25*Rf
    !!    see more datail in KBR95 Section 2.c
    rf  = r(rrmax) ! neither condition is met
    irf = rrmax
    do j = 1, prmax
       !! 1.4) calculate starting point (see KBR95 Appendix A)
       i = minval(maxloc(twd(:,j)))              !Rm:maximum tangential wind within 1200km distance
       if ( i < ira.and. irb < i ) then         !Rm does not exist within [Ra,Rb] case
          i = minval(maxloc(twd(ira:irb,j)))     ! Rm reset to relative maximum within [Ra,Rb]
          if ( minval(twd(i:irb,j)) < 0 ) then   ! anticyclonic wind exist within [Rm, Rb] case
             i = irb                             !  Rm = Rb 
          else                                   ! cyclonic within [Rm, Rb] case
             i = i                               !  Rm = Rm
          end if
       end if

       ir_strt = searchidx(r,1.1*r(i),-1) ! starting radius is set to 1.1*Rm     
       rs(j) = r(ir_strt)
       !! 2.1) calculate Rf
       n = 0
       do i = ir_strt, rrmax
          ! V < 3m/s
          if ( twd(i,j) < 3. ) then ! second condition is met
             rf(j) = r(i)
             irf(j) = i
             exit
          end if
          ! V < 6 m/s and -dV/dr < 4x10^-6 /s   !1degree ~ 100 km = 1e5 m
          if ( twd(i,j) < 6. .and. -(twd(i+1,j)-twd(i,j))/1.e5/drr < 4.e-6 ) then 
             n = n + 1
             if ( n == 2 ) then   ! first condition is met
                rf(j) = r(i)
                irf(j) = i
                exit
             end if
          end if
       end do
    end do

    !! 2.2) calculate R0
    do j = 1, prmax
       m = irf(j)
       n = min(searchidx(r,1.25*rf(j),1),rrmax)     
       if ( minval(twd(m:n,j)) < 0 ) then  !anticyconic wind exist within [Rf,1.25*Rf] case
          do i = m, n
             if (twd(i,j)<0) exit          
          end do
          r0(j) = r(i)                     !R0 set to the radius of the innermost negative V within [Rf,1.25*Rf]
       else
          r0(j) = min(1.25*rf(j),r(rrmax)) !R0 set to 1.25*Rf
       end if
    end do

    write(*,*) " determination of filter radius (degree)"
    write(*,'("  Ra, Rmax_ave, Rb, Rf_ave= ",f5.1,f5.1,f5.1,f5.1)'),ra,r_max,rb,rf_ave
    write(*,'("  NE direction Rf, R0= ",f5.1,f5.1)'),rf(4), r0(4) 
    write(*,'("  NW direction Rf, R0= ",f5.1,f5.1)'),rf(10),r0(10) 
    write(*,'("  SW direction Rf, R0= ",f5.1,f5.1)'),rf(16),r0(16)
    write(*,'("  SE direction Rf, R0= ",f5.1,f5.1)'),rf(22),r0(22) 
    n = minval(minloc(r0))
    write(*,'("  minimum      Phi,R0= ",f5.1,f5.1)'),phi(n)*180./PI,r0(n)
    n = minval(maxloc(r0))
    write(*,'("  maximum      Phi,R0= ",f5.1,f5.1)'),phi(n)*180./PI,r0(n)

    ! for debug
    if (debug) then
       open(21,file='check_domain.txt')
       write(21,*) cx,cy
       do i = 1, prmax
          write(21,*) cx+r0(i)*cos(phi(i)),cy+r0(i)*sin(phi(i))
       end do
       write(21,*) cx+r0(1)*cos(phi(1)),cy+r0(1)*sin(phi(1))
       close(21)
    end if
  end subroutine filter_domain

  subroutine remove_vortex(imax,jmax,kmax,xi,yi,dstrb,cx,cy,r0,vrtex)
    !
    ! DESCRIPTION
    !   Split disturbance field into non-hurricane disturbance field and analyzed 
    !  vortex component by application of the optimal interpolation method (Gandin 1963).
    !  In solving the lenear equation system, a non-negative least square solution is used.
    !  See more detail in KBR95 Section 2.d and Appendix B.
    !
    ! ARGUMENTS
    !  INPUT:
    !    cx,cy                 : TC center position
    !    xi(imax)              : input x-coordinate
    !    yi(jmax)              : input y-coordinate
    !    dstrb(imax,jmax,kmax) : disturbance field filtered by SMOOTHING_FILTER
    !    r0(prmax)             : filter domain radius determined by FILTER_DOMAIN
    !  OUTPUT:
    !    vrtex(imax,jmax,kmax) : analyzed vortex field
    !
    ! MODULE
    !  ip_module - grd2cyl,searchidx
    !
    ! SUBROUTINE CALLED
    !  llns - solve non-negative least square problems
    !
    ! MODULE PARAMETER
    !  integer(kind=i4b) :: prmax
    !  real(kind=sp)     :: drr,dscale
    !
    integer(kind=i4b), intent(in) :: imax,jmax,kmax
    real(kind=sp), intent(in)     :: cx,cy
    real(kind=sp), intent(in)     :: xi(imax),yi(jmax),dstrb(imax,jmax,kmax),r0(prmax)
    real(kind=sp), intent(out)    :: vrtex(imax,jmax,kmax)

    integer(kind=i4b) :: i,j,k,n,ierr
    real(kind=sp)     :: x1,x2,y1,y2,dist2,maxr0,minr0,rp,phip,r0p,dphi,novrt
    real(kind=sp)     :: r(rrmax),phi(prmax),cdstrb(rrmax,prmax,kmax),bound(prmax,kmax)
    real(kind=dp)     :: a(prmax,prmax),b(prmax),w(prmax)
    integer(kind=i4b) :: iwork(prmax)
    real(kind=dp)     :: work(prmax),rnorm

    ! cyrindrical coordinate
    do i = 1, rrmax
       r(i) = drr*i
    end do
    dphi = 2.*PI/real(prmax)
    do i = 1, prmax
       phi(i) = dphi*(i-1)
    end do

    minr0 = minval(r0)
    maxr0 = maxval(r0)

    ! 2.0) calculate correlation coefficient betweeen bouding points; mu_{i,j} in eq.(B.6) => coefficiend matrix A
    a = 0.
    do i = 1, prmax
       x1 = r0(i)*cos(phi(i))
       y1 = r0(i)*sin(phi(i))
       do j = i, prmax
          x2 = r0(j)*cos(phi(j))
          y2 = r0(j)*sin(phi(j))
          dist2 = (x1-x2)**2 + (y1-y2)**2          !dist**2
          a(i,j) = dble(exp(-dist2/dscale2))       !(dist/dscale)**2
          a(j,i) = a(i,j)
       end do
    end do

    ! 1) calculate the disturbance field value at the points on the filter radii R0
    call grd2cyl(imax,jmax,kmax,xi,yi,cx,cy,dstrb,rrmax,prmax,r,phi,cdstrb)
    do j = 1, prmax
       i = searchidx(r,r0(j),0)
       bound(j,:)  =  cdstrb(i,j,:)
    end do

    vrtex = 0.
    ! 2) calculate non-hurricane disturbance field by the optimum interpolation method
    do i = 1, imax
       x1 = xi(i) - cx
       do j = 1, jmax                    
          y1 = yi(j) - cy
          
          ! 2.0) check whether (x1,y1) is within filter domain
          rp = sqrt(x1**2+y1**2)
          if (rp > maxr0) cycle
          if (rp > minr0) then
             phip = atan2(y1,x1)
             if (phip < 0) phip = phip + 2.*PI   ! 0<=phip<2*PI
             n = searchidx(phi,phip,1)           
             if (n == prmax) then                ! phi(prmax)<=phip<phi(1)
                r0p = (r0(1)-r0(prmax))/dphi * (phip-phi(prmax)) + r0(prmax)
             else                                ! phi(n)<=phip<phi(n+1)
                r0p = (r0(n+1)-r0(n))/dphi * (phip-phi(n)) + r0(n)
             end if
             if (rp > r0p) cycle
          end if

          ! if (xi,yi) is within filter domain, the optimal interpolation is applied
          ! 2.1) correlation coefficient; mu_{p,i} in eq.(B.6) => R.H.S. vector b
          do n = 1, prmax
             x2 = r0(n)*cos(phi(n))
             y2 = r0(n)*sin(phi(n))
             dist2 = (x1-x2)**2 + (y1-y2)**2 
             b(n) = dble(exp(-dist2/dscale2))
          end do

          ! 2.2) solve linear equation system Aw=b
          !      w_{ij,n}mu_{m,n} = mu_{ij,m} (m=1,phimax)       
          call nnls(a,prmax,prmax,b,w,rnorm,work,iwork,ierr)
          if ( ierr /= 1 ) then
             print*, "error lapack routine DGELS. ierr=", ierr
             call abort
          end if

          ! 2.3) calculate non-hurricane componet of disturbance field by eq.(B.3) 
          do k = 1, kmax
             novrt = 0.
             do n = 1, prmax          
                if (w(n) < 0) then
                   print*, xi(i),yi(j),k,w(n)
                   call abort
                end if
                novrt = novrt + real(w(n))*bound(n,k)
             end do
             ! 2.4) vortex disturbance componet
             vrtex(i,j,k) = dstrb(i,j,k) - novrt
          end do
        end do
    end do
 
  end subroutine remove_vortex
  

end module vortex_module


