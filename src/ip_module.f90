module ip_module
  !
  ! DESCRIPTION
  !  this module contains interpolation routines
  !
  use type_module
  use constant_module, only : PI=>math_pi
  implicit none

  private ::
  public :: searchidx, interp2d, grd2cyl

contains

  subroutine interp2d(imaxi,jmaxi,kmax,xi,yi,datai,imaxo,jmaxo,xo,yo,datao,undef_value)
    !
    ! DESCRIPTION
    !  bilinear horizontal interpolation
    !
    ! ARGUMENTS
    !  INPUT: 
    !    xi(imaxi)               : input grid x-coordinate
    !    yi(jmaxi)               : input grid y-coordinate
    !    datai(imaxi,jmaxi,kmax) : input global grid value
    !    xo(imaxo)               : output grid x-coordinate
    !    xo(jmaxo)               : output grid y-coordinate
    !    undef_value             : undefiend value of the point out of input grid (defalt=9.999e+20)
    !  OUTPUT:
    !    datao(imaxo,jmaxo,kmax) : output global grid value
    !
    !  NOTE:
    !   output filed outside of input grid are set to 'undef_value'
    !
    integer(kind=i4b), intent(in)  :: imaxi,jmaxi,kmax,imaxo,jmaxo
    real(kind=sp),     intent(in)  :: xi(imaxi),yi(jmaxi),datai(imaxi,jmaxi,kmax)
    real(kind=sp),     intent(in)  :: xo(imaxo),yo(jmaxo)
    real(kind=sp),     intent(out) :: datao(imaxo,jmaxo,kmax)
    real(kind=sp), intent(in), optional :: undef_value
    integer(kind=i4b) :: i,j,k,ix,jy
    real(kind=sp) :: dx,dy

    if (present(undef_value)) then
       datao = undef_value
    else
       datao = 9.999e+20
    end if
    do j = 1, jmaxo          
       jy = searchidx(yi,yo(j),1)                 !yi(jy) => yo -> yi(jy+1)
       if ( jy==0 .or. jy==jmaxi+1) cycle
       dy = (yo(j)-yi(jy))/(yi(jy+1)-yi(jy))      !dy:(1-dy)
       do i = 1, imaxo
          ix = searchidx(xi,xo(i),1)              !xi(ix) => xo -> xi(ix+1)
          if ( ix==0 .or. ix==imaxi+1) cycle
          dx = (xo(i)-xi(ix))/(xi(ix+1)-xi(ix))   !dx:(1-dx)        
          ! bilinear interpolation
          do k = 1, kmax
             datao(i,j,k) = (1-dy) * ( (1-dx)*datai(ix,jy,k) + dx*datai(ix+1,jy,k) ) &
                  &       + dy * ( (1-dx)*datai(ix,jy+1,k) + dx*datai(ix+1,jy+1,k) ) 
          end do
       end do
    end do

  end subroutine interp2d
  subroutine grd2cyl(imax,jmax,kmax,xi,yi,cx,cy,grd,rmax,pmax,r,phi,cyl)
    !
    ! DESCRIPTION
    !  interpolate grid data from catersian coordinate to 
    !  cylindrical coordinate centered at (cx,cy)
    !
    ! ARGUMENTS
    !  INPUT:
    !    xi(imax)  : x value of catersian coordinate
    !    yi(jmax)  : y value of catersian coordinate
    !    r(rmax)   : radial value of cyrindrical coordinate
    !    phi(pmax) : azimuthal value of cyrindrical coordinate (RADIAN)
    !    grd(imax,jmax,kmax) : catersian grid value
    !    cx,cy     : center position in catersian coordinate
    !  OUTPUT    
    !    cyl(rmax,pmax,kmax) : cyrindrical grid value
    !
    integer(kind=i4b), intent(in)  :: imax,jmax,kmax,rmax,pmax
    real(kind=sp),     intent(in)  :: cx,cy,xi(imax),yi(jmax),r(rmax),phi(pmax)
    real(kind=sp),     intent(in)  :: grd(imax,jmax,kmax)
    real(kind=sp),     intent(out) :: cyl(rmax,pmax,kmax)
    
    integer(kind=i4b) :: i,j,k,ix,jy
    real(kind=sp) :: x,y,dx,dy,ex,ey,xx(imax),yy(jmax)
    
    ! coordinate relative to (cx,cy)
    xx = xi - cx 
    yy = yi - cy
    do j = 1, pmax
       ex = cos(phi(j))
       ey = sin(phi(j))        
       do i = 1, rmax      
          x  = r(i)*ex
          ix = searchidx(xx,x,1)
          dx = (x-xx(ix))/(xx(ix+1)-xx(ix))      !dx:(1-dx)

          y  = r(i)*ey
          jy = searchidx(yy,y,1)
          dy = (y-yy(jy))/(yy(jy+1)-yy(jy))      !dy:(1-dy)
          
          if (ix<1 .or. ix>imax .or. jy<1 .or. jy>jmax) then
             print*, "the request domain is out of the bound"
             print*, "request (x,y) is ",x+cx,y+cy
             print*, "input domain (xmin,xmax,ymin,ymax) is",xi(1),xi(imax),yi(1),yi(jmax)
             call abort
          end if
          do k = 1, kmax             
             cyl(i,j,k) = (1-dy) * ( (1-dx)*grd(ix,jy,k) + dx*grd(ix+1,jy,k) ) &
                  &     + dy * ( (1-dx)*grd(ix,jy+1,k) + dx*grd(ix+1,jy+1,k) )       
          end do
       end do
    end do
    
  end subroutine grd2cyl
  function searchidx(arr, val, opt) result(idx)
    !
    ! DESCRIPTION
    !  search for the index of array 
    !
    ! ARGUMENTS
    !  INPUT:
    !    arr : 1-dimensional array, which must be monotinically increasing or decreasing
    !    val : value to search for 
    !    opt : option flag
    !      opt>0 - arr(idx) <= var < arr(idx+1) or arr(idx) >= var > arr(idx+1)
    !      opt<0 - arr(idx+1) < var <= arr(idx) or arr(idx+1) > var >= arr(idx)
    !      opt=0 - return the first index nearest to the value 'val'
    !  OUTPUT:
    !    idx : index of the array
    !
    real(kind=sp), intent(in)      :: arr(:), val
    integer(kind=i4b), intent(in)  :: opt
    integer(kind=i4b) :: idx
    integer(kind=i4b) :: imax 
    
    imax = size(arr)
    idx = minval(minloc(abs(arr-val)))
    if (opt>0) then        
       if (arr(1) < arr(imax)) then !increasing case          
          if (arr(idx) > val) idx = idx - 1          
       else                         !decreasing case
          if (arr(idx) < val) idx = idx - 1          
       end if
    else if (opt<0) then   
       if (arr(1) < arr(imax)) then !increasing case          
          if (arr(idx) < val) idx = idx + 1          
       else                         !decreasing case
          if (arr(idx) > val) idx = idx + 1          
       end if
    end if

  end function searchidx

end module ip_module


