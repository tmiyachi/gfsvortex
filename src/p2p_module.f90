module p2p_module
  !
  ! DESCRIPTION
  !  this module contains vertical interpolation module
  !
  use type_module
  use ip_module
  use calcmet_module
  implicit none

  private
  public :: p2p, p2p_extrapolate_T, p2p_q2rh
  
contains  

  subroutine p2p(imax,kmaxi,datai,pi,ps,kmaxo,po,datao,constant_value)
    !
    ! DESCRIPTION 
    !   vertically interpolate input filed from pressure level to another pressue level
    !  by 4-point lagrange interpolation method.
    !   
    ! NOTE
    !   1) The value of layers adjacent to the input top and bottom layer are calculated by
    !     linear interpolation. Abobe the highest input level and under the lowest input level
    !     the field is kept constant and is equatl to the value of the highest input level
    !     and the lowest input level respectively. 
    !
    !                  constant (use top level field)
    !     pi(kmaxi)    ---------------------------------------
    !                  linear (use 2 levels)
    !     pi(kmaxi-1)  ---------------------------------------
    !                  lagrange interpolation (use 4 levels)
    !     pi(2)        ---------------------------------------
    !                  linear (use 2 levels)
    !     pi(1)        --------------------------------------- 
    !                  constant_value
    !   
    !   2) If input pressure level is under the surface, the field at the level is not used
    !     for the interpolation.
    !
    !                  lagrange interpolation (use 4 levels)
    !     pi(kps+1)    -----------------------------------------
    !                  linear (use 2 levels)
    !     pi(kps)      -----------------------------------------
    !     ps           -----------------------------------------
    !                  constant_value
    !
    ! ARGUMENTS
    !  INPUT:
    !   datai(imax,kmaxi) : input field to interpolate
    !   pi(imax,kmaxi)    : input field pressure level
    !   ps(imax)          : input surface pressure field
    !   po(imax,kmaxo)    : output filed pressure level
    !   constant_value    : constant value of the field under the bottom level (optional)
    !  OUTPUT:  
    !   datao(imax,kmaxo) : interpolated field
    !    
    ! MODULE
    !  ip_module - searchidx
    !
    integer(kind=i4b), intent(in)  :: imax,kmaxi,kmaxo
    real(kind=sp),     intent(in)  :: datai(imax,kmaxi),ps(imax)
    real(kind=sp),     intent(in)  :: pi(imax,kmaxi)
    real(kind=sp),     intent(in)  :: po(imax,kmaxo)
    real(kind=sp),     intent(out) :: datao(imax,kmaxo)
    real(kind=sp),intent(in),optional :: constant_value

    integer(kind=i4b) :: i,k,ki,kps
    real(kind=sp) :: lnpo
    real(kind=sp) :: lnpi(imax,kmaxi)
    real(kind=sp) :: p0,p1,p2,p3,p4,q1,q2,q3,q4,a

    lnpi = -log(pi)

    !$omp parallel private(k,ki,kps,lnpo,a,p0,p1,p2,p3,p4,q1,q2,q3,q4)
    !$omp do
    do i = 1, imax
       kps  = searchidx(lnpi(i,:),-log(ps(i)),-1) !lnpi(kps-1)<=-ln(ps)<lnpi(kps)
       do k = 1, kmaxo
          lnpo = -log(po(i,k))
          ki = searchidx(lnpi(i,:),lnpo,1)        !lnpi(ki)<=lnpo<lonpi(ki+1)
          if ( ki == 0 .or.  ki < kps ) then            !!out of bottom level
             if (present(constant_value)) then
                datao(i,k) = constant_value
             else
                datao(i,k) = datai(i,ki+1)    !constant
             end if
          else if ( ki == kmaxi) then                   !!out of top level
             datao(i,k) = datai(i,kmaxi)      !constant
          else if ( ki == kps .or. ki+1 == kmaxi ) then !!bottom or top level
             a = (datai(i,ki+1)-datai(i,ki)) / (lnpi(i,ki+1)-lnpi(i,ki))
             datao(i,k) = a*(lnpo-lnpi(i,ki)) + datai(i,ki) !linear interpolation
          else                                          !!p1<p2<po<p3<p4
             p0 = lnpo
             p1 = lnpi(i,ki-1)
             p2 = lnpi(i,ki)
             p3 = lnpi(i,ki+1)
             p4 = lnpi(i,ki+2)

             q1 = (p0-p2)*(p0-p3)*(p0-p4)/(p1-p2)/(p1-p3)/(p1-p4)
             q2 = (p0-p1)*(p0-p3)*(p0-p4)/(p2-p1)/(p2-p3)/(p2-p4)
             q3 = (p0-p1)*(p0-p2)*(p0-p4)/(p3-p1)/(p3-p2)/(p3-p4)
             q4 = (p0-p1)*(p0-p2)*(p0-p3)/(p4-p1)/(p4-p2)/(p4-p3)
             
             datao(i,k) = datai(i,ki-1)*q1 + datai(i,ki)*q2 +  &
                  &       datai(i,ki+1)*q3 + datai(i,ki+2)*q4     !cubic lagrangian

          end if
       end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine p2p

  subroutine p2p_extrapolate_T(imax,kmaxi,tempi,pi,psi,kmaxo,po,pso,hs,tempo)
    ! DESCRIPTION 
    !   vertically interpolate input filed from pressure level to another pressue level
    !  by 4-point lagrange interpolation method for "temperature field"
    !   
    ! NOTE
    !     The layers adjacent to the input top and bottom layer are calculated using
    !     linear interpolation. Above the highest input level the temperature is kept constant
    !     and equal to the value of the highest input level. Between the lowest input level and 
    !     the input surface the temperature is interpolated using EXTRAPOLATE_T.
    !
    !                  constant (equal to the value of highest level)
    !     pi(kmaxi)    ----------------------------------------------
    !                  linear (use 2 levels)
    !     pi(kmaxi-1)  ----------------------------------------------
    !                  lagrange interpolation (use 4 levels)
    !     pi(2)        ----------------------------------------------
    !                  linear (use 2 levels)
    !     pi(1)        ----------------------------------------------
    !                  linearly extrapolated by EXTRAPOLATE_T
    !     ps           ----------------------------------------------
    !                  linearly extrapolated by EXTRAPOLATE_T
    !
    ! ARGUMENTS
    !  INPUT:
    !   tempi(imax,kmaxi) : input temperature field to interpolate [K]
    !   pi(imax,kmaxi)    : input field pressure level             [Pa]
    !   psi(imax)         : input surface pressure field           [Pa]
    !   hs(imax)          : output surface height                  [m]
    !   po(imax,kmaxo)    : output filed pressure level            [Pa]
    !   pso(imax)         : output pressure surface field          [Pa]
    !  OUTPUT:  
    !   tempo(imax,kmaxo) : interpolated temperature field
    !
    ! MODULE
    !  ip_module - searchidx
    !
    integer(kind=i4b), intent(in)  :: imax,kmaxi,kmaxo
    real(kind=sp),     intent(in)  :: tempi(imax,kmaxi),psi(imax),pso(imax),hs(imax)
    real(kind=sp),     intent(in)  :: pi(imax,kmaxi),po(imax,kmaxo)
    real(kind=sp),     intent(out) :: tempo(imax,kmaxo)

    integer(kind=i4b) :: i,k,ki,kps
    real(kind=sp) :: lnpo,Ts
    real(kind=sp) :: lnpi(imax,kmaxi)
    real(kind=sp) :: p0,p1,p2,p3,p4,q1,q2,q3,q4,a

    lnpi = -log(pi)

    !$omp parallel private(k,ki,kps,lnpo,a,p0,p1,p2,p3,p4,q1,q2,q3,q4,Ts)
    !$omp do
    do i = 1, imax
       kps  = searchidx(lnpi(i,:),-log(psi(i)),-1) !lnpi(kps-1)<=-ln(psi)<lnpi(kps)
       Ts = extrapolate_Ts(tempi(i,kps),pi(i,kps),pso(i)) !calculate surface temperature of output field
       do k = 1, kmaxo
          lnpo = -log(po(i,k))
          ki = searchidx(lnpi(i,:),lnpo,1)        !lnpi(ki)<=lnpo<lonpi(ki+1)
          if (po(i,k) > pso(i)) then              !!under surface
             !extrapolate from surface             
             tempo(i,k) = extrapolate_T(po(i,k),pso(i),hs(i),Ts) 
          else if ( ki == 0 .or.  ki < kps ) then       !!out of bottom level          
             !extrapolate from surface             
             tempo(i,k) = extrapolate_T(po(i,k),pso(i),hs(i),Ts)
          else if ( ki == kmaxi) then                   !!out of top level
             tempo(i,k) = tempi(i,kmaxi)      !constant
          else if ( ki == kps .or. ki+1 == kmaxi ) then !!bottom or top level
             a = (tempi(i,ki+1)-tempi(i,ki)) / (lnpi(i,ki+1)-lnpi(i,ki))
             tempo(i,k) = a*(lnpo-lnpi(i,ki)) + tempi(i,ki) !linear interpolation
          else                                          !!p1<p2<po<p3<p4
             p0 = lnpo
             p1 = lnpi(i,ki-1)
             p2 = lnpi(i,ki)
             p3 = lnpi(i,ki+1)
             p4 = lnpi(i,ki+2)

             q1 = (p0-p2)*(p0-p3)*(p0-p4)/(p1-p2)/(p1-p3)/(p1-p4)
             q2 = (p0-p1)*(p0-p3)*(p0-p4)/(p2-p1)/(p2-p3)/(p2-p4)
             q3 = (p0-p1)*(p0-p2)*(p0-p4)/(p3-p1)/(p3-p2)/(p3-p4)
             q4 = (p0-p1)*(p0-p2)*(p0-p3)/(p4-p1)/(p4-p2)/(p4-p3)
             
             tempo(i,k) = tempi(i,ki-1)*q1 + tempi(i,ki)*q2 +  &
                  &       tempi(i,ki+1)*q3 + tempi(i,ki+2)*q4     !cubic lagrangian

          end if
       end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine p2p_extrapolate_T
  
  subroutine p2p_q2rh(imax,kmaxi,qi,tempi,pi,psi,kmaxo,po,pso,hs,qo)
    !
    ! DESCRIPTION
    !  vertically interpolate specific humidity field to output pressire level
    !  Under the lowest input level and the input surface, the specific humidity
    !  is extrapolated to be kept constant value of the relative humidity.
    !
    ! ARGUMENT
    !  INPUT:
    !   qi(imax,kmaxi)    : input specfic humidity field
    !   tempi(imax,kmaxi) : input temperature field
    !   pi(imax,kmaxi)    : input field pressure level
    !   psi(imax)         : input surface pressure field
    !   po(imax,kmaxo)    : output filed pressure level  
    !   pso(imax)         : output surface pressure field
    !   hs(imax,kmaxo)    : geometric height
    !  OUTPUT:  
    !   qo(imax,kmaxo)    : interpolated specific humidity field
    !
    integer(kind=i4b), intent(in)  :: imax,kmaxi,kmaxo
    real(kind=sp),     intent(in)  :: qi(imax,kmaxi),tempi(imax,kmaxi)
    real(kind=sp),     intent(in)  :: pi(imax,kmaxi),psi(imax)
    real(kind=sp),     intent(in)  :: po(imax,kmaxo),pso(imax),hs(imax)
    real(kind=sp),     intent(out) :: qo(imax,kmaxo)

    integer(kind=i4b) :: i, k
    real(kind=sp)     :: rhi(imax,kmaxi), tempo(imax,kmaxi)

    do k = 1, kmaxi
       do i = 1, imax
          rhi(i,k) = calcmet_q2rh(tempi(i,k),qi(i,k),pi(i,k))
       end do
    end do

    call p2p(imax,kmaxi,rhi,pi,psi,kmaxo,po,qo)
    call p2p_extrapolate_T(imax,kmaxi,tempi,pi,psi,kmaxo,po,pso,hs,tempo)
    qo = max(min(qo,100.),0.)

    do k = 1, kmaxi
       do i = 1, imax
          qo(i,k) = calcmet_rh2q(tempo(i,k),qo(i,k),po(i,k))
       end do
    end do

  end subroutine p2p_q2rh


end module p2p_module
